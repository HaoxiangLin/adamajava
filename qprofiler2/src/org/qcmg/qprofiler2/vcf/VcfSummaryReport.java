package org.qcmg.qprofiler2.vcf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;

import org.apache.commons.lang3.StringUtils;
import org.qcmg.common.model.ProfileType;
import org.qcmg.common.util.QprofilerXmlUtils;
import org.qcmg.common.vcf.VcfFormatFieldRecord;
import org.qcmg.common.vcf.VcfRecord;
import org.qcmg.common.vcf.header.VcfHeader;
import org.qcmg.qprofiler2.report.SummaryReport;
import org.qcmg.qprofiler2.summarise.SampleSummary;
import org.qcmg.qprofiler2.util.XmlUtils;
import org.w3c.dom.Element;


public class VcfSummaryReport  extends SummaryReport {
	
//	public static final String NodeHeader = "vcfHeader";
//	public static final String NodeHeaderMeta = "MetaInformation" ;
//	public static final String NodeHeaderMetaLine = "MetaInformationLine";
//	public static final String NodeHeaderStructured = "StructuredMetaInformation";
//	public static final String NodeHeaderStructuredType = "StructuredMetaInformationType";
//	public static final String NodeHeaderStructuredLine = "StructuredMetaInformationLine";
//	public static final String NodeHeaderFinal =  "HeaderLine";		
//	public static final String NodeCategory = "ReportingCategory";
//	public static final String id = "id";
	
	public static final String Sample = "sample";	
//	static final String Seperator = ",:,";
	static final String Seperator = ":";
	private final VcfHeader vcfHeader;	
	private final String[] sampleIds; 
 
	//it allows the format field value eg. --formart FT=PASS, then it seperate value to PASS the others
	private final String[] formatCategories;
//	Map< String, SampleSummary > summaries = new HashMap<>();
	//Map< sample, Map<joined_cat,SampleSummary> > summaries
	Map< String, Map<String,SampleSummary> > summaries = new HashMap<>();
	Map< String, AtomicLong > counts = new HashMap<>();
	
	public VcfSummaryReport(VcfHeader header, String[] formats){	     
		this.vcfHeader = header;
		this.formatCategories = formats; 
		this.sampleIds = header.getSampleId(); 	
	}
	
	public void toXml(Element parent) {
		logger.info("preparing output...");
		Element parentElement = init(parent, ProfileType.VCF);	
		logger.info("outputing vcf header to xml...");
		XmlUtils.vcfHeaderToXml(parentElement, vcfHeader);
		logger.info("outputing sample information to xml...");
		summaryToXml( parentElement  );		
	}
	

	void parseRecord( VcfRecord  vcf ) {
		updateRecordsParsed();
		
		List<String> formats = vcf.getFormatFields();
		if(sampleIds == null || formats.size() != sampleIds.length + 1)
			logger.warn("missing/redundant sample column exists in vcf record: " + vcf.toSimpleString());
		
		//for each sample column
		for(int i = 1; i < formats.size(); i ++){
		//	String key = sampleIds[i-1]; 
			VcfFormatFieldRecord re = new VcfFormatFieldRecord(formats.get(0), formats.get(i));
			 
			List<String> cates = new ArrayList<>();
			for(String cate : formatCategories ){
				//new
				int pos = cate.indexOf("="); 
				if(pos > 0){
					String formatKey = cate.substring(0, pos).trim();
					String formatValue = cate.substring(pos+1).trim();
					cates.add(  re.getField(formatKey) == null ? null :
							  re.getField(formatKey).equalsIgnoreCase(formatValue) ? formatValue : "Other"  );						
				} else {					 
					cates.add(   re.getField(cate) );
					 
				}
			}	 
			Map<String, SampleSummary> map =  summaries.computeIfAbsent( sampleIds[i-1], (k) -> new HashMap<String, SampleSummary>() );
			map.computeIfAbsent( StringUtils.join( cates, Seperator ), (k) -> new SampleSummary() ).parseRecord( vcf, i ) ;
		}				
	}

	/**
	 * modifying now
	 * @param parent
	 */
	void summaryToXml(Element parent){			
		//get list of types eg. FT:INF:CONF
		List<String>  formatsTypes = new ArrayList<>();
		for(int i = 0; i < formatCategories.length; i ++) {
			int pos = formatCategories[i].indexOf("=");
			formatsTypes.add( pos > 0 ? formatCategories[i].substring(0, pos) : formatCategories[i] );			
		}	
		
		Element summaryElement =  QprofilerXmlUtils.createSubElement(parent,  ProfileType.VCF.getReportName()+"Metrics" );		
		for( String sample : summaries.keySet() ) {	
			Element ele =  QprofilerXmlUtils.createSubElement( summaryElement, Sample);
			ele.setAttribute(XmlUtils.Sname, sample);
						
			//Element ele = getSampleElement( summaryElement, sample );
			for(String cates : summaries.get(sample).keySet() )			
				if( formatsTypes.isEmpty() )
					summaries.get(sample).get(cates).toXML( ele, null, null );
				else
					summaries.get(sample).get(cates).toXML( ele, String.join(Seperator , formatsTypes), cates );				
		}		
	}

//	private Element getSampleElement(Element parent, String sampleId) {
//		 
//		List<Element> Esamples =  QprofilerXmlUtils.getOffspringElementByTagName(parent, Sample);
//		for(Element ele : Esamples)  
//			if( ele.getAttribute(XmlUtils.Sname ).equals(sampleId))
//				return ele; 
//		 
//		Element ele =  QprofilerXmlUtils.createSubElement( parent, Sample);
//		ele.setAttribute(XmlUtils.Sname, sampleId);
//		
//		return ele;
//	}
	
}
