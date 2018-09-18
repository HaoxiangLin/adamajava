package org.qcmg.qprofiler2.summarise;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

import org.qcmg.common.model.QCMGAtomicLongArray;
import org.qcmg.common.model.SubsitutionEnum;
import org.qcmg.common.string.StringUtils;
import org.qcmg.common.util.Constants;
import org.qcmg.common.util.IndelUtils;
import org.qcmg.common.util.IndelUtils.SVTYPE;
import org.qcmg.common.util.Qprofiler1XmlUtils;
import org.qcmg.common.vcf.VcfFormatFieldRecord;
import org.qcmg.common.vcf.VcfRecord;
import org.qcmg.common.vcf.header.VcfHeaderUtils;
import org.qcmg.qprofiler2.util.XmlUtils;
import org.w3c.dom.Element;

public class SampleSummary {
	public final static String variantType = "VariantType";
	public final static String genotypes = "Genotypes";
	public final static String genotype = "Genotype";
	public final static String substitutions = "Substitutions";	
	public final static String substitution = "Substitution";
	public final static String report = "Report";	
	public final static String VAF = "VariantAlleleFrequency";
	public final static String binTally = "BinTally";
	public final static String tiTvRatio = "TiTvRatio";
	public final static String transitions ="Transitions";
	public final static String transversions = "Transversions";
		
	//allele frequency value range. eg. altBinSize = 10, each bin contain counts for variants which mutation allele percent is [0, 0.1]
	public final static int altBinSize = 20; 
		
	Set< String > gts = new HashSet<>(); //store possible genotyp 0/1, 0/0, ./1 ...
	Map< String, AtomicLong > summary = new HashMap<>();
	Map<String, QCMGAtomicLongArray> summaryAD = new HashMap<>();
	
	AtomicLong counts = new AtomicLong();
	QCMGAtomicLongArray sampleTrans = new QCMGAtomicLongArray(  SubsitutionEnum.values().length );	

	private void increment(String key){		
		summary.computeIfAbsent(key, v -> new AtomicLong()).incrementAndGet();
	}
	
	/**
	 * Updates the supplied map with 
	 *  
	 * NOT SIDE-EFFECT FREE
	 * 
	 * @param type
	 * @param gt
	 * @param ad
	 * @param dp
	 * @param map
	 */
	public static void incrementGTAD(SVTYPE type,String gt, String ad, String dp, Map<String, QCMGAtomicLongArray> map){	
		
		if( ad == null || ad.contains(".") ||  gt == null || gt.contains(".")  ||  gt.equals("0/0") || gt.equals("0|0") ) 
				return;		
		
		int commaIndex = ad.indexOf(Constants.COMMA);
		int vaf = 0;
		
		/*
		 * vaf needs to equal the sum of all the numbers in the AD field, apart from the first number (which is the reference count)
		 */
		while (commaIndex > -1) {
			int nextConmmaIndex = ad.indexOf(Constants.COMMA, commaIndex + 1);
			if (nextConmmaIndex > -1) {
				vaf += Integer.parseInt(ad.substring(commaIndex + 1, nextConmmaIndex));
				commaIndex = nextConmmaIndex;
			} else {
				vaf += Integer.parseInt(ad.substring(commaIndex + 1));
				break;
			}
		}
		 
		int rate = (int) ( 0.5 + (double) ( vaf * altBinSize ) / Integer.parseInt(dp) );
		map.computeIfAbsent(type.name() , (k) -> new QCMGAtomicLongArray(altBinSize + 1)).increment( rate );	
	}
		

	public void parseRecord( VcfRecord  vcf, int formateOrder) {
		VcfFormatFieldRecord format = new VcfFormatFieldRecord(vcf.getFormatFields().get(0), vcf.getFormatFields().get(formateOrder));
		counts.incrementAndGet(); //total number
		SVTYPE type = IndelUtils.getVariantType(vcf.getRef(), vcf.getAlt());
		boolean isdbSNP = !StringUtils.isNullOrEmptyOrMissingData(vcf.getId());
		
		String gt = format.getField(VcfHeaderUtils.FORMAT_GENOTYPE); //GT
		gts.add(gt);
				
		increment( type + ""); //count svtype
		increment( type + gt);	//count genotyp	
						
		//variant allel frequence VAF	 
		incrementGTAD(type , gt, format.getField(VcfHeaderUtils.FORMAT_ALLELIC_DEPTHS), format.getField(VcfHeaderUtils.FORMAT_READ_DEPTH), summaryAD);				
 								
		if(isdbSNP) increment(  type.name() + "dbSNP");	//dbsnp					
		if(type.equals(SVTYPE.SNP)){ 
			//get Ti Tv counts based on GT value			
			String salt = vcf.getAlt().replace(",", "");
			String sgt = gt.replace("|", "").replace("/", "").replace(".", "").replace("0", "");				
			List<SubsitutionEnum> transTypes = new ArrayList<>();
			new HashSet<Character> (sgt.chars().mapToObj(e->(char)e).collect(Collectors.toList())).forEach(
				c ->  transTypes.add( SubsitutionEnum.getTrans( vcf.getRef().charAt(0), salt.charAt( c-'1' ) )) 
			);
						
			for(SubsitutionEnum transType : transTypes){
			//	String mark = transType.isTranstion() ? "Ti" : transType.isTransversion() ? "Tv" : "Other";	
			//	increment( type.name() + mark  );
				increment( type.name() + transType.name() );					
			}				
		}		
	}	 			
		
	public void toXML(Element parent, String formats, String values ){
		Element reportE = Qprofiler1XmlUtils.createSubElement(parent, report);
		if(formats != null) {
			reportE.setAttribute("foramts", formats);
			reportE.setAttribute("values", values);
		}
		List<String> orderedGTs = new ArrayList<>(gts);
		orderedGTs.sort(null);
		for(SVTYPE type : SVTYPE.values()){	
			//only output non zero value
			AtomicLong totalAL = summary.get( type.name());
			if( null == totalAL) continue;
			
			Element svtypeE = Qprofiler1XmlUtils.createSubElement(reportE,  variantType);
			svtypeE.setAttribute("type", type.toVariantType() );
			svtypeE.setAttribute( XmlUtils.count, totalAL.toString());
			String key =  type.name() + "dbSNP";
			svtypeE.setAttribute("inDBSNP", summary .containsKey(key)? summary.get(key).get()+"" : "0" ); 
						
			Map<String, AtomicLong> gtvalues = new HashMap<>();
			for(String gt : orderedGTs) {
				AtomicLong gtv =  summary.get(type.name() + gt);
				if(gtv != null) gtvalues.put(gt, gtv);				
			}

			XmlUtils.outputMap(svtypeE, "variant counts", "genotype distribution based on Format column GT field value", "genotype4"+type.toVariantType(), "genotyp(GT field)", gtvalues);

			QCMGAtomicLongArray array = summaryAD.get(type.name());
			if(array != null) {
				Map<String, AtomicLong> altFreq = new HashMap<>();					
				for( int i = 0; i < altBinSize; i ++  ){
					long count =  array.get(i) ;
					if(count <= 0) continue;
					
					String  bin = String.format("%.2f,%.2f", (double) i / altBinSize, (double) (i+1) / altBinSize  );
					altFreq.put(bin, new AtomicLong(count) );
				}		
							
				XmlUtils.outputMap(svtypeE, "variant counts", "variant distribution based on allele frequency value", "altFreq4"+type.toVariantType(), "variantAltFrequency rate range", altFreq);
			}  
			
			//titv
			if(type.equals(SVTYPE.SNP)){
				Map<String, AtomicLong> tiFreq = new HashMap<>();
				Map<String, AtomicLong> tvFreq = new HashMap<>();
				for(SubsitutionEnum tran: SubsitutionEnum.values()) 
					if(tran.isTranstion() &&  summary.get(type.name()+tran.name()) != null) 
						tiFreq.put(tran.toString(), summary.get(type.name()+tran.name()));
					else if( tran.isTransversion() &&  summary.get(type.name()+tran.name()) != null)  					 
						tvFreq.put(tran.toString(), summary.get(type.name()+tran.name()));			
				XmlUtils.outputMap(svtypeE, "variant counts", "SNP transtion substitution distribution", "transtion" , "substitution change", tiFreq);
				XmlUtils.outputMap(svtypeE, "variant counts", "SNP transversion substitution distribution", "transversion" , "substitution change", tvFreq);
			}			
		

		}
	}
}
