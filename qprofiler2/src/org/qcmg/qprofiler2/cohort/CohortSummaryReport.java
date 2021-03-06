package org.qcmg.qprofiler2.cohort;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

import org.qcmg.common.util.IndelUtils.SVTYPE;
import org.qcmg.common.util.QprofilerXmlUtils;
import org.qcmg.qprofiler2.report.SummaryReport;
import org.qcmg.qprofiler2.summarise.SampleSummary;
import org.qcmg.qprofiler2.util.XmlUtils;
import org.qcmg.qprofiler2.vcf.VcfSummaryReport;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

public class CohortSummaryReport extends SummaryReport {
	public static String outputSeperate = "\t";	
	public static String headerline = String.join( outputSeperate, new String[]{ "File","Sample" ,"ReportingCategory","VariantType","VariantCount","DbSnpProportion","TiTvRatio"});
	
	final String sampleId; 
	final String file;	
	List<Category> categories = new ArrayList<>();
	
	private long sum_ti = 0; 
	private long sum_tv = 0;
	private long sum_db = 0;
	private long sum_count = 0;
	
	CohortSummaryReport(File fxml, Element sampleNode) throws IOException{
		this.file = fxml.getCanonicalPath();
		this.sampleId = sampleNode.getAttribute(XmlUtils.Sid);
		
		for( Element ele : QprofilerXmlUtils.getChildElementByTagName( sampleNode, SampleSummary.report ) ) {
			Category cat = new Category( ele.getAttribute(SampleSummary.values), ele ); 
			categories.add(cat );
			
			//summary
			sum_ti += cat.ti;
			sum_tv += cat.tv;
			for(int db : cat.dbSnpCounts.values()) sum_db += db;
			for(int count : cat.variantsCounts.values()) sum_count += count;			
		}
				
	}	
	
	List<String> outputCounts(  ){		 
		List<String> output = new ArrayList<>(); 
		for(Category cat: categories)
			for(String str : cat.output())
				output.add(file + outputSeperate + sampleId + outputSeperate + str);
		
		return output;
	}	
	
	/**
	 * summary counts will output to log file
	 * @return an array of total counts of variantCounts , dbsnp , ti and tv
	 */
	public long[] getCountSum(){
		return new long[]{ sum_count, sum_db, sum_ti, sum_tv  };		
	}
	
	class Category{
		final String category; //eg. FORMAT:FT=PASS;FORMAT:INF=.
		HashMap<String, Integer> variantsCounts = new HashMap<>();
		HashMap<String, Integer> dbSnpCounts = new HashMap<>();
		
		//for snv only, there is only one svn for each report unde same category
		final String titvRate; 
		int ti = 0, tv = 0;
						
		Category(String name, Element report){
					 
			this.category = (name == null || name.isEmpty()) ? "-" : name; 
			String titv = "-" ;
			 
			//for(Element ele :QprofilerXmlUtils.getOffspringElementByTagName(report, SampleSummary.variantType)){
			for(Element ele :QprofilerXmlUtils.getOffspringElementByTagName(report, XmlUtils.metricsEle)){
				//record counts and dbsnp for all type variants
				String type = ele.getAttribute(XmlUtils.Sname);	
				int count = Integer.parseInt(ele.getAttribute("count"));
				variantsCounts.put(type, count );
				
				Element e1 = QprofilerXmlUtils.getChildElementByTagName(ele, XmlUtils.Svalue).stream().filter(e -> e.getAttribute(XmlUtils.Sname).equals(SampleSummary.dbSNP)).findFirst().get();				
				int db = Integer.parseInt( e1.getTextContent() );
				dbSnpCounts.put(type, db);
				
				if(type.equals(SVTYPE.SNP.toVariantType())){
					e1 = QprofilerXmlUtils.getChildElementByTagName(ele, XmlUtils.Svalue).stream().filter(e -> e.getAttribute(XmlUtils.Sname).equals(SampleSummary.tiTvRatio)).findFirst().get();				
					titv = e1.getTextContent();	
					
					//ti
					Optional<Element> streams = QprofilerXmlUtils.getChildElementByTagName(ele, XmlUtils.variableGroupEle).stream().filter(e -> e.getAttribute(XmlUtils.Sname).equals(SampleSummary.transitions)).findFirst();				
					if( streams.isPresent()   ) {
						List<Integer> sums = new ArrayList<>();						
						QprofilerXmlUtils.getChildElementByTagName(streams.get(), XmlUtils.Stally).stream()
							.forEach( e ->   sums.add(  Integer.parseInt(e.getAttribute(XmlUtils.Scount)) ) );
						ti = sums.stream().mapToInt(i -> i.intValue()).sum();
					} 
					
					//tv
					streams = QprofilerXmlUtils.getChildElementByTagName(ele, XmlUtils.variableGroupEle).stream().filter(e -> e.getAttribute(XmlUtils.Sname).equals(SampleSummary.transversions)).findFirst() ;				
					if( streams.isPresent() ) {
						List<Integer> sums = new ArrayList<>();
						QprofilerXmlUtils.getChildElementByTagName(streams.get(), XmlUtils.Stally).stream()
							.forEach( e ->   sums.add(  Integer.parseInt(e.getAttribute(XmlUtils.Scount)) ) );
						tv = sums.stream().mapToInt(i -> i.intValue()).sum();
					}
					
				}
			}
			 
			this.titvRate = titv;				 
		}
		
		/**
		 * order: "ReportingCategory","VariantType","VariantCount","DbSnpProportion","TiTvRatio"
		 * @return a list of counts for each variant type
		 */
		List<String> output(){
			List<String> output = new ArrayList<>(); 
			for(String type :  variantsCounts.keySet() ){	
				StringBuilder sb = new StringBuilder(category);
				sb.append(outputSeperate).append(type)
				  .append(outputSeperate).append( variantsCounts.get(type) )
				  .append(outputSeperate).append( String.format( "%.3f",  (double) dbSnpCounts.get(type) / variantsCounts.get(type) ))
				  .append(outputSeperate).append( type.equals(SVTYPE.SNP.toVariantType()) ? titvRate : "-" );
				
				output.add( sb.toString() );
			}
			return output;
		}		
	}

	
	@Override @Deprecated
	public void toXml(Element parent) {}

}
