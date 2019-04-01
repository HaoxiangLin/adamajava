package org.qcmg.qprofiler2.summarise;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import org.qcmg.common.log.QLogger;
import org.qcmg.common.log.QLoggerFactory;
import org.qcmg.common.model.QCMGAtomicLongArray;
import org.qcmg.common.util.QprofilerXmlUtils;
import org.qcmg.qprofiler2.summarise.PairSummary.Pair;
import org.qcmg.qprofiler2.util.XmlUtils;
import org.w3c.dom.Element;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadGroupSummary {
//	protected QLogger logger = QLoggerFactory.getLogger(getClass());
	private long errNo = 0 ;
	public final static int errReadLimit  = 10;	
	//xml node name 	
	public final static String node_readgroup = "readGroup";		
	public final static String node_softClip = "softClippedBases";
	public final static String node_trim = "trimmedBases";	
	public final static String node_hardClip = "hardClippedBases";
	public final static String node_readLength = "readLength" ; 
	public final static String node_overlap = "overlappedBases";	
	public final static String node_duplicate = "duplicateReads";
	public final static String node_secondary = "secondary";
	public final static String node_supplementary = "supplementary"; 	
	public final static String node_unmapped = "unmappedReads";
	public final static String node_nonCanonicalPair = "nonCanonicalPairs";
	public final static String node_failedVendorQuality = "failedVendorQuality";
	public final static String modal_isize = "modalSize";
		
	public final static String smin= "min";	
	public final static String smax = "max";
	public final static String smean = "mean"; 
	public final static String smode =  "mode"; 
	public final static String smedian = "median" ; 
	public final static String stdDev = "stdDev" ; 
	
	public final static String sreadCount = "readCount"; 
	public final static String slostBase = "basesLost"; 
	public final static String sbasePercent = "basePercent"; 
	public static final String sbaseCount = "baseCount"; //"countedBase"
	
 			
	//softclips, hardclips, read length; 	
	QCMGAtomicLongArray softClip = new QCMGAtomicLongArray(128);
	QCMGAtomicLongArray hardClip = new QCMGAtomicLongArray(128);
	QCMGAtomicLongArray readLength = new QCMGAtomicLongArray(128);		
//	QCMGAtomicLongArray isize = new QCMGAtomicLongArray(PairSummary.middleTlenValue);	 //store count bwt [0, 1499]
	private final ConcurrentMap<String, AtomicLong> cigarValuesCount = new ConcurrentHashMap<String, AtomicLong>();
	//must be concurrent set for multi threads
//	Set<PairSummary> pairCategory = Collections.newSetFromMap(new ConcurrentHashMap<PairSummary, Boolean>());
	private final ConcurrentMap<Integer, PairSummary> pairCategory = new ConcurrentHashMap<>();

	AtomicInteger max_isize = new AtomicInteger(); 	
	//bad reads inforamtion
	AtomicLong duplicate = new AtomicLong();
	AtomicLong secondary  = new AtomicLong();
	AtomicLong supplementary  = new AtomicLong();
	AtomicLong unmapped  = new AtomicLong();
	AtomicLong unpaired  = new AtomicLong();
	
//	AtomicLong nonCanonical  = new AtomicLong();
	AtomicLong failedVendorQuality  = new AtomicLong();
	AtomicLong inputReadCounts  = new AtomicLong();
	
	
	//for combined readgroups, since max read length maybe different
	//below value will be reset on preSummary() inside readsSummary2Xml()
	private long trimedBase = -1;
	private long maxBases = -1; 
	private long duplicateBase = -1;
	private long unmappedBase = -1;
	private long nonCanonicalBase = -1;

	private int maxReadLength = 0;
	private long noOfCountedReads = -1;
	private long trimedRead = -1;
	
		
	private final String readGroupId; 		
	public ReadGroupSummary(String rgId){
		this.readGroupId = rgId; 
	}
	
	public String getReadGroupId(){
		return readGroupId; 
	}	
	
	public Long getMaxBases(){	
		return this.maxBases; 
	}

	/**
	 * 
	 * @return read base of duplicated, unmapped, non-canonical paired 
	 */
	public long getDuplicateBase() { 
		return readGroupId .equals( QprofilerXmlUtils.All_READGROUP ) ? 
			this.duplicateBase : this.duplicate.get() * this.maxReadLength;
	}
	
	public long getUnmappedBase() {
		return readGroupId .equals( QprofilerXmlUtils.All_READGROUP ) ? 
			this.unmappedBase :  this.unmapped.get() * this.maxReadLength;
	}
	
	public long getnonCanonicalBase() {
		return readGroupId .equals( QprofilerXmlUtils.All_READGROUP ) ? 
			this.nonCanonicalBase : getnonCanonicalReadsCount() * this.maxReadLength;
	}
	/**
	 * 
	 * @return the sum of reads marked as not proper pair
	 */
	public long getnonCanonicalReadsCount() {
		long sum = 0;
		for(PairSummary p : pairCategory.values()) {
			if(p.isProperPair == false) {
				sum += p.getFirstOfPairCounts();
				sum += p.getSecondOfPairCounts();
			}
		}
		return sum; 		
	}
	
	public long getTrimmedBase() {
		return trimedBase;
	}
	
	public long getTrimmedRead() {
		return trimedRead;
	}
	

		
	/**
	 * classify record belongs (duplicate...) and count the number. If record is parsed, then count the hard/soft clip bases, pair mapping overlap bases and pairs information
	 * @param record
	 * @return true if record parsed; otherwise return false since record duplicate, supplementary, secondary, failedVendorQuality, unmapped or nonCanonical. 
	 */
	public boolean parseRecord( final SAMRecord record ){
				
		//record input reads number
		inputReadCounts.incrementAndGet();  
				
		//find discard reads and return false
		if ( record.getSupplementaryAlignmentFlag()) {
			supplementary.incrementAndGet();
			return false;
		}else if( record.isSecondaryAlignment() ) {
			secondary.incrementAndGet();
			return false;
		}else if(record.getReadFailsVendorQualityCheckFlag()) {
			failedVendorQuality.incrementAndGet();
			return false;
		} 
		
		//parseing cigar
		//cigar string from reads including duplicateReads, nonCanonicalPairs and unmappedReads but excluding discardedReads (failed, secondary and supplementary).
		parseCigar(record.getCigar());		
					
		//find the max read length	
		if(record.getReadLength() > this.maxReadLength) {
			this.maxReadLength = record.getReadLength();			
		}
		
		//find mapped badly reads and return false	
		if(record.getDuplicateReadFlag()){
			duplicate.incrementAndGet();
			return false;
		}


		//check pair orientaiton, tLen, mate
		if(record.getReadPairedFlag()) {
			if(record.getReadUnmappedFlag() && record.getMateUnmappedFlag()) {				 
				unmapped.incrementAndGet();
				return false;
			}
			PairSummary.Pair pairType = PairSummary.getPairType(record);
			boolean isProper = record.getProperPairFlag();

			int key = isProper? pairType.id * 2 :  pairType.id;
			if(!pairCategory.containsKey(key))
				pairCategory.put(key, new PairSummary( pairType, isProper));				
			PairSummary oldP = pairCategory.get(key);//PairSummary.computeIfAbsent(pairCategory, record);			 
			oldP.parse(record);					
		} else {
			unpaired.incrementAndGet();
		    if(record.getReadUnmappedFlag() ){
				unmapped.incrementAndGet();
				return false;
			} 
		}
					
		
		//parse clips
		 int lHard = 0, lSoft = 0;
		 for (CigarElement ce : record.getCigar().getCigarElements()) {
			 if (ce.getOperator().equals(CigarOperator.HARD_CLIP)) {
				 lHard += ce.getLength();
			 } else if (ce.getOperator().equals(CigarOperator.SOFT_CLIP)) {
				 lSoft += ce.getLength();
			 }
		 }
		 
		hardClip.increment(lHard);
		softClip.increment(lSoft);
		//filtered reads and discard reads won't count the length
		int lwithHard = record.getReadLength()+lHard;
		readLength.increment(lwithHard);
		if(lwithHard > this.maxReadLength)
			this.maxReadLength = lwithHard;
		 				
		return true;  
	}
	//debug
	private void parseCigar(Cigar cigar) {
		if (null == cigar)  return;
		
		for (CigarElement ce : cigar.getCigarElements()) {
			CigarOperator operator = ce.getOperator();
			if ( ! CigarOperator.M.equals(operator)) {
				String key = "" + ce.getLength() + operator;
				cigarValuesCount.computeIfAbsent(key, k -> new AtomicLong(0)).incrementAndGet();					
			}			 
		}	 
	}
	
	public ConcurrentMap<String, AtomicLong> getCigarCount() {		 
		return cigarValuesCount;
	}
	
	public QCMGAtomicLongArray getOverlapCount(){
		QCMGAtomicLongArray overlapBase = new QCMGAtomicLongArray(PairSummary.segmentSize);	
		
		for(PairSummary p : pairCategory.values()) {
			if( !p.isProperPair) continue; 
			for(int i = 0; i < PairSummary.segmentSize; i ++) {			 	
				overlapBase.increment(i, p.getoverlapCounts().get(i) );					 
			}			
		}						
		return overlapBase;
	}	

			
	public long getDiscardreads() {
		return supplementary.get() + failedVendorQuality.get() + secondary.get();
	}
		
	/**
	 * 
	 * @return number of reads excluds discarded one. 
	 */
	public long getCountedReads() {				
		long totalRead = 0  ;
		for (int i = 1 ; i < readLength.length() ; i++)	{		 
			totalRead += readLength.get(i);
		}
		return totalRead + duplicate.get() + unmapped.get() + getnonCanonicalReadsCount() ;		
	}
	
	public int getAveReadLength(){
		long totalgoodRead = 0 , totalgoodBase = 0;
		for (int i = 1 ; i < readLength.length() ; i++){ 			 
			totalgoodBase += i * readLength.get(i);
			totalgoodRead += readLength.get(i);
		}		
		return ( totalgoodRead == 0 )? 0 :(int) (totalgoodBase / totalgoodRead);		
	}
		
	/**
	 * check all globle value and assign the sumamry value
	 * eg.  
	 * 	private long trimedBase = 0; 	
	 */
	public void preSummary(long max , long trimB, long trimR, long baseDup, long baseUnmap, long baseNonCan) {
		//maxReadLength is continuelly updated during parseRecord
				
		//number of reads excluds discarded one.
		this.noOfCountedReads = getCountedReads();
		
		//add value to All_READGROUP manually since the readLength are different for each readGroup
		if( readGroupId.equals(QprofilerXmlUtils.All_READGROUP) ) {	 			
			this.maxBases = max;
			this.trimedBase = trimB;
			this.trimedRead = trimR;
			this.duplicateBase = baseDup;
			this.unmappedBase = baseUnmap;
			this.nonCanonicalBase = baseNonCan;			
			return;
		}
		
		this.maxBases = this.noOfCountedReads * this.maxReadLength;
		
		//suppose original reads with same length for same read group
		//the counts should be same to the sumOf QCMGAtomicLongArray from parseTrim(...)
		long totalBase = 0;		
		long no = 0;
		for (int i = 1 ; i < this.maxReadLength ; i++){ 			 
			totalBase += i * readLength.get(i);
			no += readLength.get(i);
		}				
		this.trimedRead = no; 
		//here we omit full length reads which length is maxReadLength since there is no trimming
		this.trimedBase = no * this.maxReadLength - totalBase;
		
	}
	private void checkPrepare() throws Exception {
		if( !readGroupId.equals(QprofilerXmlUtils.All_READGROUP) )
			preSummary(0, 0,0,0,0,0);
		else if (trimedBase < 0 || trimedRead < 0 || maxBases < 0 ||   noOfCountedReads < 0)  
			throw new Exception("counts have null value, please call method::preSummary(...) before process!");		
	}
	
	public void readSummary2Xml(Element parent ) throws Exception { 	
		
		 checkPrepare();						
		//add to xml RG_Counts
		Element rgElement = XmlUtils.createMetricsNode(parent,"reads",  inputReadCounts.get());					
		//add discarded read Stats to readgroup summary		
		Element ele = XmlUtils.createGroupNode(rgElement, QprofilerXmlUtils.discardReads );
		XmlUtils.outputValueNode(ele, "supplementaryAlignmentCount", supplementary.get());
		XmlUtils.outputValueNode(ele, "secondaryAlignmentCount", secondary.get());
		XmlUtils.outputValueNode(ele, "failedVendorQualityCount",failedVendorQuality.get()  );
					
		//long noOfRecords = getCountedReads( );
		badReadStats( rgElement, node_duplicate, duplicate.get(), getDuplicateBase()   );
		badReadStats( rgElement, node_unmapped, unmapped.get(), getUnmappedBase() );		
		badReadStats( rgElement, node_nonCanonicalPair, getnonCanonicalReadsCount(),  getnonCanonicalBase());
				
		long lostBase = getDuplicateBase() + getUnmappedBase()+ getnonCanonicalBase();
		
		if( !readGroupId.equals(QprofilerXmlUtils.All_READGROUP) ){
			lostBase += lostBaseStats( rgElement, node_trim, parseTrim(readLength ), this.maxBases);				
		}else{
			lostBase += this.trimedBase;
			badReadStats( rgElement, node_trim, this.trimedRead, this.trimedBase   );
		}
					
		lostBase += lostBaseStats( rgElement, node_softClip, softClip, this.maxBases );
		lostBase += lostBaseStats( rgElement, node_hardClip, hardClip, this.maxBases );
		lostBase += lostBaseStats( rgElement, node_overlap, getOverlapCount(), this.maxBases );
		XmlUtils.addCommentChild(( Element )rgElement.getLastChild(), "Only count overlapped bases on strand for pairs which have a positive TLEN value.");
		
		//create node for overall
		Element overallEle = XmlUtils.createGroupNode(rgElement, QprofilerXmlUtils.overall );
		XmlUtils.outputValueNode( overallEle, "readMaxLength", this.maxReadLength  );
		XmlUtils.outputValueNode( overallEle, "readAveLength",getAveReadLength());
		//readCount
		String comment = sreadCount + ": includes duplicateReads, nonCanonicalPairs and unmappedReads but excludes discardedReads (failed, secondary and supplementary).";
		overallEle.appendChild( overallEle.getOwnerDocument().createComment(comment) );				
		XmlUtils.outputValueNode( overallEle, sreadCount,  this.noOfCountedReads );
		
		XmlUtils.outputValueNode( overallEle, "readUnpairedCount",  unpaired.get() );
		
		
		//baseCount
		comment = sbaseCount +  (readGroupId.equals(QprofilerXmlUtils.All_READGROUP)? 
				": the sum of " + sbaseCount + " from all read group" 	: 	": " + sreadCount + " * readMaxLength");
		overallEle.appendChild( overallEle.getOwnerDocument().createComment(comment) );
		XmlUtils.outputValueNode( overallEle, sbaseCount, this.maxBases  );		
		
		//baseLost
		comment = slostBase +  (readGroupId.equals(QprofilerXmlUtils.All_READGROUP)? 
				": the sum of " + slostBase + " from all read group" 	: 	": readMaxLength * (duplicateReads + nonCanonicalPairs + unmappedReads) + trimmedBases + softClippedBases + hardClippedBases + overlappedBases");
		overallEle.appendChild( overallEle.getOwnerDocument().createComment(comment) );		
		XmlUtils.outputValueNode( overallEle, slostBase,  lostBase);	
		
				
		//add overall information to current readgroup element	
		comment = String.format("%s: %s / %s", QprofilerXmlUtils.lostPercent, slostBase, sbaseCount);
		overallEle.appendChild( overallEle.getOwnerDocument().createComment(comment) );
		double lostPercent =  this.maxBases == 0? 0: 100 * (double) lostBase / this.maxBases ;			
		XmlUtils.outputValueNode( overallEle, QprofilerXmlUtils.lostPercent, lostPercent );			
	}
	 	 
	public void pairSummary2Xml( Element parent ) { 
		//add to xml RG_Counts
		Element ele =  XmlUtils.createMetricsNode( parent, "properPairs", null );
		long sum = 0;
		for(PairSummary p : pairCategory.values()) {
			if( p.isProperPair) {
				p.toSummaryXml(ele);	
				sum += p.getPairCounts();
			}			
		}
		ele.setAttribute( XmlUtils.Scount, sum+"");  
		
		ele =  XmlUtils.createMetricsNode( parent, "notProperPairs", null );
		sum = 0;
		for(PairSummary p : pairCategory.values()) {
			if(! p.isProperPair) {
				p.toSummaryXml(ele);
				sum += p.getPairCounts();
			}			
		}	
		ele.setAttribute( XmlUtils.Scount, sum+"");  
		
			 
//		//output tlen refer to method lostBaseStats(...)
//		long sum = 0, no = 0;
//		for(int i = 1; i < PairSummary.middleTlenValue; i ++ ){
//			sum += isize.get(i) * i;
//			no += isize.get(i);
//		}	
//		
//		int mean = (no == 0) ? 0: (int) (sum / no);		
//		int mode = 0; //mode is the number of read which length is most popular
//		long highest = 0;
//		double sd = 0;
//		for (int i = 1 ; i < PairSummary.middleTlenValue ; i++) {	
//			//the number of isize.get(i) pairs have same isize value 
//			sd += Math.pow(( i - mean), 2) * isize.get(i) / no;
//			if(isize.get(i) > highest){ highest = isize.get(i); mode = i;  } 
//		}
//		double standardDeviation = Math.sqrt(sd);
//		
//		//create node for overall, tLen
//		Element overallEle = XmlUtils.createGroupNode(ele, QprofilerXmlUtils.overall );
//		XmlUtils.outputValueNode(overallEle, "maxTlen", this.max_isize.get());
//		XmlUtils.outputValueNode(overallEle, "pairCount", ( f5f3.getPairCounts() +  f3f5.getPairCounts() + inward.getPairCounts() + outward.getPairCounts() ) );	
//
//		overallEle.appendChild(overallEle.getOwnerDocument().createComment("below stats is based on tlen value < " + PairSummary.middleTlenValue));
//		XmlUtils.outputValueNode(overallEle, "meanUnderTlen5000", mean);
//		XmlUtils.outputValueNode(overallEle, "modeUnderTlen5000", mode);	
//		XmlUtils.outputValueNode(overallEle, "pairCountUnderTlen5000", no);		
//		XmlUtils.outputValueNode( overallEle, "stdDevUnderTlen5000", (int)standardDeviation);			

	}
	
	public void pairTlen2Xml( Element parent ) {		
		for(PairSummary p : pairCategory.values()) {
			if( p.isProperPair) {
				p.toTlenXml(parent);				
			}			
		}
	}
	
	private void badReadStats(Element parent, String nodeName, long reads, long badBase ){
		Element ele = XmlUtils.createGroupNode(parent, nodeName);				
		XmlUtils.outputValueNode(ele, sreadCount, reads);	
		XmlUtils.outputValueNode(ele, slostBase, badBase);		
		double percentage = 100 * (double) badBase / this.maxBases ;
		XmlUtils.outputValueNode(ele, "basePercent",  percentage);		 	 
	}	
			
	private long lostBaseStats(Element parent, String nodeName, QCMGAtomicLongArray array, long maxBases ){
		long arrayLength = null != array ? array.length() : 0;
		
		long bases = 0,counts = 0;		
		for (int i = 1 ; i < arrayLength ; i++){
			if(array.get(i) <= 0) continue;
			counts += array.get(i);
			bases += i * array.get(i);
		}		
		int mean = (counts == 0) ? 0: (int) (bases / counts);	
				
		 // to avoid aray.get(0) >= 0 since 1(counts)/2== 0(counts/2) == 0
		long medium = 0;  	
		for (int i = 1 ; i < arrayLength; i++) {
			if(( medium += array.get(i)) > counts/2 ){ medium = i;  break; }
		}
		int min = 0; //find the smallest non-zero value;
		for(int i = 1; i < arrayLength; i ++) {
			if(array.get(i) > 0){ 
				min  = i; break; 
			}
		}
		
		int max = 0; //find the biggest non-zero value;
		for(int i = (int) (arrayLength -1); i > 0; i--) {
			if(array.get(i) > 0){ 
				max = i; break;  
			}
		}
		
		int mode = 0; //mode is the number of read which length is most popular
		long highest = 0;
		for (int i = 1 ; i < arrayLength ; i++) { 					
			if(array.get(i) > highest){
				highest = array.get(i);
				mode = i; 
			}  	
		}
		Element ele = XmlUtils.createGroupNode(parent, nodeName);	
		XmlUtils.outputValueNode(ele, smin, min);
		XmlUtils.outputValueNode(ele, smax, max);
		XmlUtils.outputValueNode(ele, smean, mean);
		XmlUtils.outputValueNode(ele, smode, mode);	
		XmlUtils.outputValueNode(ele, smedian, medium);
		XmlUtils.outputValueNode(ele, sreadCount, counts);
		XmlUtils.outputValueNode(ele, slostBase, bases);
		
		//deal with boundary value, missing reads
		double percentage = (maxBases == 0)? 0: 100 * (double) bases /  maxBases ;				
		XmlUtils.outputValueNode(ele, QprofilerXmlUtils.lostPercent,  percentage );	
									
		return bases; 
	}
	
	private QCMGAtomicLongArray parseTrim( QCMGAtomicLongArray readLengthArray ) {
		QCMGAtomicLongArray array = new QCMGAtomicLongArray(  maxReadLength+1 ); 
		for(int i = 1; i < maxReadLength; i ++) {
			array.increment( maxReadLength - i, readLengthArray.get(i) );
		}
		return array;
	}



}