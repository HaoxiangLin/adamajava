/**
 * © Copyright The University of Queensland 2010-2014.  This code is released under the terms outlined in the included LICENSE file.
 */
package org.qcmg.qsv;

import java.io.File;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicInteger;

import org.qcmg.gff3.GFF3FileReader;
import org.qcmg.gff3.GFF3Record;
import org.qcmg.qsv.discordantpair.DiscordantPairCluster;
import org.qcmg.qsv.discordantpair.PairGroup;
import org.qcmg.qsv.report.DCCReport;
import org.qcmg.qsv.report.QSVClusterReport;

/**
 * Class to write QSV cluster details to output files
 * @author felicity
 *
 */
public class QSVClusterWriter {
   
    private QSVParameters tumorParameters;
    private QSVParameters normalParameters;
    private AtomicInteger somaticCount = new AtomicInteger(0);
	private AtomicInteger germlineCount = new AtomicInteger(0);
    private AtomicInteger normalGermlineCount = new AtomicInteger(0);
	private boolean isQCMG;
	private String analysisId;
	private boolean isSingleSided;
	private boolean twoFileMode;
	private int minInsertSize;
	private String validationPlatform;
	private List<String> gffFiles;
	private Map<String, List<GFF3Record>> gffMap;
    
    public QSVClusterWriter(QSVParameters tumor, QSVParameters normal, boolean isQCMG, String analysisId, boolean singleSided, boolean twoFileMode, int minInsertSize, String validationPlatform, List<String> gffFiles) throws Exception {
    	 this.tumorParameters = tumor;
         this.normalParameters = normal;
         this.isQCMG = isQCMG;
         this.analysisId = analysisId;
         this.isSingleSided = singleSided;
         this.twoFileMode = twoFileMode;
         this.minInsertSize = minInsertSize;
         this.validationPlatform = validationPlatform;
         this.gffFiles = gffFiles;         
         this.gffMap = parseGFFFiles();
         
	}
	
	private Map<String, List<GFF3Record>> parseGFFFiles() throws Exception {
		Map<String, List<GFF3Record>> gffMap = new HashMap<String, List<GFF3Record>>();
		for (String file: gffFiles) {
			GFF3FileReader reader = new GFF3FileReader(new File(file));
			Iterator<GFF3Record> it = reader.getRecordIterator();
			while (it.hasNext()) {
				GFF3Record g3 = it.next();
				if (gffMap.containsKey(g3.getSeqId())) {
					gffMap.get(g3.getSeqId()).add(g3);
				} else {
					List<GFF3Record> list = new ArrayList<GFF3Record>();
					list.add(g3);
					gffMap.put(g3.getSeqId(), list);					
				}
			}
			reader.close();
		}
		
		return gffMap;		
	}

	/**
	 * Write cluster records
	 * @param clusterRecords
	 * @param isTumour
	 * @throws Exception
	 */
	public synchronized void writeQSVClusterRecords(Map<PairGroup, Map<String, List<DiscordantPairCluster>>> clusterRecords, boolean isTumour) throws Exception {
		List<QSVCluster> svRecords = new ArrayList<QSVCluster>();
    	if (clusterRecords != null) {
    		//any discordant pair records that aren't already converted to QSV cluster records
	        for (Entry<PairGroup, Map<String, List<DiscordantPairCluster>>> entry: clusterRecords.entrySet()) {
	            for (Map.Entry<String, List<DiscordantPairCluster>> mutationTypeEntry : entry.getValue().entrySet()) {   
                    List<DiscordantPairCluster> currentClusters = mutationTypeEntry.getValue();               
                    
                    for (DiscordantPairCluster r: currentClusters) {                    	
                    	QSVCluster record = new QSVCluster(r, false, tumorParameters.getSampleId());
                    	svRecords.add(record);
                    }                    
	            }
	        }       
    	}
    	
    	if (isTumour) {
    		writeTumourSVRecords(svRecords);
    	} else {
    		writeNormalSVRecords(svRecords);
    	}
	}

	/**
	 * Write the records for the tumour/test sample to file 
	 * @param svRecords
	 * @throws Exception
	 */
	public synchronized void writeTumourSVRecords(List<QSVCluster> svRecords) throws Exception {
    	String base = tumorParameters.getResultsDir();
    	Date analysisDate = tumorParameters.getAnalysisDate();
    	String sampleId = tumorParameters.getSampleId();
    	
    	List<QSVCluster> somaticRecords = new ArrayList<QSVCluster>();
    	List<QSVCluster> germlineRecords = new ArrayList<QSVCluster>();
    	
    	for (QSVCluster record: svRecords) {    		
    		String id = "";
    		
    		 if (record.printRecord(isSingleSided)) {
    			 
    			 if (record.passesMinInsertSize(minInsertSize)) {
    				
    				if (twoFileMode) {
			    		if (record.isGermline()) {
			    			if (isQCMG) {
			    				id += "stgv_" + germlineCount.incrementAndGet();
			    			} else {
			    				id += "gm_" + germlineCount.incrementAndGet();
			    			}
			    		} else {
			    			if (isQCMG) {
			    				id += "stsm_" + somaticCount.incrementAndGet();
			    			} else {
			    				id += "sm_" + somaticCount.incrementAndGet();
			    			}
			    		}
    				} else {
    					
    					if (isQCMG) {
    						id += "stgv_" + somaticCount.incrementAndGet();
    					} else {
    						id += "sv_" + somaticCount.incrementAndGet();
    					}
    				}
    				
    				record.checkReferenceFlank(tumorParameters.getReference(), tumorParameters.getChromosomes());
    				record.checkGFF(gffMap);
    				
		    		
		    		String svId = id;
		    		record.setIdParameters(svId, analysisId, sampleId, analysisDate);
		    		
		    		if (record.isGermline()) {		    			
		    			germlineRecords.add(record);
		    		} else {		    			
		    			somaticRecords.add(record);    			
		    		}
    			 }
    		 }
	   	}
    	
    	writeReports(base, "somatic", somaticRecords, analysisDate, sampleId);
    	writeReports(base, "germline", germlineRecords, analysisDate, sampleId);
	}	

	/**
	 * Write the records for the normal/control sample to file 
	 * @param svRecords
	 * @throws Exception
	 */
	public synchronized void writeNormalSVRecords(List<QSVCluster> svRecords) throws Exception {
		String type = "normal-germline";
        String base = normalParameters.getResultsDir();
        Date analysisDate = normalParameters.getAnalysisDate();
    	String sampleId = normalParameters.getSampleId();
    	  
    	if (this.isQCMG) {
		    	for (QSVCluster record: svRecords) {
		    		
		    		 if (record.printRecord(isSingleSided)) {
		    			 if (record.passesMinInsertSize(minInsertSize)) {
		    				 
		    				 String svId = "stnmgv_" + normalGermlineCount.incrementAndGet();
		    				 record.setIdParameters(svId,analysisId, sampleId, analysisDate);
		    			 }
		    		 }
		    	}
		    	
		    	writeReports(base, type, svRecords, analysisDate, sampleId);				
    	}
	}	

	/**
	 * Write reports for a QSV Cluster. Reports can be dcc, tab delimited, qprimer or verbose
	 * @param base
	 * @param type
	 * @param records
	 * @param analysisDate
	 * @param sampleId
	 * @throws Exception
	 */
	public synchronized void writeReports(String base, String type, List<QSVCluster> records, Date analysisDate, String sampleId) throws Exception {		
		String outType = type + ".";
		if (!twoFileMode) {
			outType = "";
		}
		//QCMG specific reports: dcc, qprimer, soft clip
		if (this.isQCMG) {
			if (twoFileMode) {	
				new DCCReport(new File(base + "." + outType + "dcc"), records, tumorParameters.getFindType(), normalParameters.getFindType(), isSingleSided, validationPlatform);
			} else {
				new DCCReport(new File(base + "." + outType + "dcc"), records, tumorParameters.getFindType(), null, isSingleSided, validationPlatform);
			}
			if (type.equals("somatic")) {
				new QSVClusterReport(base, outType, "qprimer", "qprimer", records, isSingleSided);
				new QSVClusterReport(base, outType, "softclip.txt", "softclip", records, isSingleSided);
			}
		}
		
		//Tab delimited
		if (twoFileMode) {			
			new QSVClusterReport(base, outType, "sv.txt", "tab", records,  tumorParameters.getFindType(), normalParameters.getFindType(), isSingleSided);			
		 } else {
			new QSVClusterReport(base, outType, "sv.txt", "tab", records, tumorParameters.getFindType(), null, isSingleSided);			
		 }
		
		//Verbose report
		for (QSVCluster r : records) {		
			if (r.printRecord(isSingleSided)) {	
				 if (!twoFileMode) {					
					new QSVClusterReport(base, outType, "records", "verbose", r, tumorParameters.getFindType(), null, isQCMG, isSingleSided);
				 } else {
 					new QSVClusterReport(base, outType, "records", "verbose", r, tumorParameters.getFindType(), normalParameters.getFindType(), isQCMG, isSingleSided);
				 }
 			}			
		}			
	}
	
    public AtomicInteger getSomaticCount() {
		return somaticCount;
	}

	public void setSomaticCount(AtomicInteger somaticCount) {
		this.somaticCount = somaticCount;
	}

	public AtomicInteger getGermlineCount() {
		return germlineCount;
	}

	public void setGermlineCount(AtomicInteger germlineCount) {
		this.germlineCount = germlineCount;
	}

	public AtomicInteger getNormalGermlineCount() {
		return normalGermlineCount;
	}

	public void setNormalGermlineCount(AtomicInteger normalGermlineCount) {
		this.normalGermlineCount = normalGermlineCount;
	} 
}
