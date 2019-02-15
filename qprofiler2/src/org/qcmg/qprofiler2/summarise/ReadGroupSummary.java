package org.qcmg.qprofiler2.summarise;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import org.qcmg.common.model.QCMGAtomicLongArray;
import org.qcmg.common.util.QprofilerXmlUtils;
import org.qcmg.picard.util.PairedRecordUtils;
import org.qcmg.qprofiler2.util.XmlUtils;
import org.w3c.dom.Element;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadGroupSummary {
	//xml node name 	
	public final static String node_readgroup = "readGroup";
	public final static String node_f5f3 = "f5f3Pair";
	public final static String node_f3f5 = "f3f5Pair";
	public final static String node_inward = "inwardPair";
	public final static String node_outward = "outwardPair";		
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
	//debug
	//for sequenceMetrics::"paris"/variableGroup::tlen
	public final static String smean = "meanUnderTlen5000"; 
	public final static String smode =  "modeUnderTlen5000"; 
	public final static String smedian = "medianUnderTlen5000" ; 
	public final static String stdDev = "stdDevUnderTlen5000" ; 
	public final static String pairCount = "pairCountUnderTlen5000" ; 
	
	public final static String sreadCount = "readCount"; 
	public final static String slostBase = "basesLost"; 
	public final static String sbasePercent = "basePercent"; 
	public static final String sbaseCount = "baseCount"; //"countedBase"
	
	//fixed value
	public final static int bigTlenValue = 10000;
	public final static int smallTlenValue = 1500;
	public final static int middleTlenValue = 5000;
	public final static int rangeGap = 100;
 			
	//softclips, hardclips, read length; 	
	QCMGAtomicLongArray softClip = new QCMGAtomicLongArray(128);
	QCMGAtomicLongArray hardClip = new QCMGAtomicLongArray(128);
	QCMGAtomicLongArray readLength = new QCMGAtomicLongArray(128);
	QCMGAtomicLongArray overlapBase = new QCMGAtomicLongArray(128);	
		 
	//QCMGAtomicLongArray.get(arrayTlenLimit) for tlen=[bigTlenValue, ~)
	QCMGAtomicLongArray isize1 = new QCMGAtomicLongArray(bigTlenValue + 1);	//previous method
	QCMGAtomicLongArray isize = new QCMGAtomicLongArray(middleTlenValue);	 //store count bwt [0, 1499]
	QCMGAtomicLongArray isizeRange = new QCMGAtomicLongArray( (bigTlenValue/rangeGap) + 1);	//store count bwt [0/100, 10000/100]
	AtomicInteger max_isize = new AtomicInteger(); 
	
	//bad reads inforamtion
	AtomicLong duplicate = new AtomicLong();
	AtomicLong secondary  = new AtomicLong();
	AtomicLong supplementary  = new AtomicLong();
	AtomicLong unmapped  = new AtomicLong();
	AtomicLong nonCanonical  = new AtomicLong();
	AtomicLong failedVendorQuality  = new AtomicLong();
	AtomicLong inputReadCounts  = new AtomicLong();
	
	//pairing information	
	AtomicLong pairNum  = new AtomicLong();
	PairedRead f3f5 = new PairedRead( node_f3f5 );
	PairedRead f5f3 = new PairedRead( node_f5f3 );
	PairedRead inward = new PairedRead( node_inward );
	PairedRead outward = new PairedRead( node_outward );	
	
	//pair in different reference
	AtomicLong diffRef  = new AtomicLong();	
	AtomicLong mateUnmapped  = new AtomicLong();	
	AtomicLong diffRef_flag_p  = new AtomicLong();	
	
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
			this.nonCanonicalBase : this.nonCanonical.get() * this.maxReadLength;
	}
	
	public long getTrimmedBase() {
		return trimedBase;
	}
	
	public long getTrimmedRead() {
		return trimedRead;
	}
	
 	private class PairedRead {
 		private final String name;		
 		PairedRead(String name){this.name = name;}
 		
		AtomicLong overlap = new AtomicLong();
		AtomicLong near = new AtomicLong();
		AtomicLong far = new AtomicLong();
		AtomicLong bigTlen  = new AtomicLong();	
		
		AtomicLong recordSum  = new AtomicLong();
		AtomicLong recordSum_flag_p = new AtomicLong();
		
		void parse(SAMRecord record, Integer overlapBase ){		
			recordSum.getAndIncrement();
			
			if(record.getProperPairFlag())	recordSum_flag_p.getAndIncrement();
			
			if( overlapBase > 0 ){
				overlap.incrementAndGet();	
			}
			if(record.getInferredInsertSize() >= bigTlenValue){
				bigTlen.incrementAndGet();
			}							
			if( record.getInferredInsertSize() >= smallTlenValue && record.getInferredInsertSize() < bigTlenValue ){
				far.incrementAndGet();		
			}
			if( record.getInferredInsertSize() < smallTlenValue && overlapBase <= 0 ){
				near.incrementAndGet();		
			}			
		}		
		
		/**
		 * 	  <f5f3Pair tlenUnder1500="3528" TlenOver10000="102075" tlenBetween1500And10000="22027" overlapping="241" count="127871"/>
	  <f3f5Pair tlenUnder1500="12784" TlenOver10000="105049" tlenBetween1500And10000="21797" overlapping="738" count="140368"/>
	  <inwardPair tlenUnder1500="134706192" TlenOver10000="198419" tlenBetween1500And10000="148179" overlapping="6988807" count="142041597"/>
	  <outwardPair tlenUnder1500="100769" TlenOver10000="205890" tlenBetween1500And10000="111947" overlapping="322649" count="741255"/>
		 * @param parent
		 */
		
		void toXml(Element parent  ){			
			Element stats = XmlUtils.createGroupNode(parent, name);
			XmlUtils.outputValueNode(stats, "overlappedPairs", overlap.get());
			//stats.appendChild(stats.getOwnerDocument().createComment("below counts excluding overlapping pairs"));
			XmlUtils.outputValueNode(stats, "tlenUnder1500Pairs", near.get() );		 
			XmlUtils.outputValueNode(stats, "tlenOver10000Pairs", bigTlen.get()  );
			XmlUtils.outputValueNode(stats, "tlenBetween1500And10000Pairs",far.get() );
			XmlUtils.outputValueNode(stats, "pairCount", recordSum.get()  );
			
		}
		
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
			
		//find the max read length	
		if(record.getReadLength() > this.maxReadLength) {
			this.maxReadLength = record.getReadLength();			
		}
		//find mapped badly reads and return false	
		if(record.getReadUnmappedFlag() ){
			unmapped.incrementAndGet();
			return false;
		}else if(record.getDuplicateReadFlag()){
			duplicate.incrementAndGet();
			return false;
		}else if(! PairedRecordUtils.isCanonical( record ) ){
			nonCanonical.incrementAndGet();
			parsePairing( record , null); 	
			return false;
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
		 		
		//parse overlap
		int overlap = PairedRecordUtils.getOverlapBase(record);
		if(overlap > 0) overlapBase.increment(overlap);
		parsePairing( record , overlap); 				
		return true; 			
 
	}
	
	public QCMGAtomicLongArray getISizeCount(){
		return isize;
	}
	public QCMGAtomicLongArray getISizeRangeCount(){
		return isizeRange;
	}
	
	/**
	 * count the first of pair number and record iSize(Tlen);
 	 * if( mate unmapped) record it and then return; 
 	 * if(mate map to different ref && first of pair) record it and then return; 
 	 * counts pair direction and distance; Here we only select reads with Tlen>0 from pair to avoiding double counts
	 * @param record: a mapped and paired read and not duplicate, not supplementary, secondary or failed
	 * @param overlapBase
	 */
	public void parsePairing( SAMRecord record, Integer overlapBase ){
		//skip non-paired reads
		if( !record.getReadPairedFlag() )  return;  
		
		//normally bam reads are mapped, if the mate is missing, we still count it to pair but no detailed pair information
		if( record.getMateUnmappedFlag() ){
			mateUnmapped.incrementAndGet(); 
			pairNum.incrementAndGet();	
			return; 
		}
		
		//pair from different reference, only look at first pair to avoid double counts
		if( !record.getReferenceName().equals( record.getMateReferenceName()) && record.getFirstOfPairFlag()){
			diffRef.incrementAndGet();	
			pairNum.incrementAndGet();			
			return; 
		}
		
		//if tlen unvaliable count first of pair and then discard	
		//pairNum >= sum(f3f5, f5f3,inward,outward, diffRef, mateUnmapped)	
		if( record.getInferredInsertSize() == 0 && record.getFirstOfPairFlag()) { 
			isize.increment(0);
			pairNum.incrementAndGet();	
			return;
		}
								
		//discard remains reads with tlen minus or unvaliable 	
		if( record.getInferredInsertSize() <= 0 ) { return; }
				
		//waringing: 
		//both pair with tlen >0 it is hardly happen, bad value but picard accept it
		//proper pair with tlen > 0 disregard first or second of pair		
		pairNum.incrementAndGet();	
		//another reason is impossible to caculate overlapBase since getOverlapBase( record) return 0 if tlen<=0		
		if(overlapBase == null)	 //non-canonical pair
			overlapBase = PairedRecordUtils.getOverlapBase( record);		
		if( PairedRecordUtils.isF5toF3(record)) f5f3.parse( record, overlapBase ); 		
		else if( PairedRecordUtils.isF3toF5(record)) f3f5.parse( record, overlapBase );	 
		else if( PairedRecordUtils.isOutward(record)) outward.parse( record, overlapBase );		
		else if( PairedRecordUtils.isInward(record)) inward.parse( record, overlapBase );
				
		
		//record tlen
		int tLen =  record.getInferredInsertSize();	
		//record the maxmum tlen value
		if( record.getInferredInsertSize() > max_isize.get() ) max_isize.getAndSet( tLen );		
		if(tLen > bigTlenValue )  tLen = bigTlenValue+1;	//can't handle too big tlen
		//record region for all pair  
		isizeRange.increment( (tLen/rangeGap));			
		//only record popular TLEN
		if(tLen < middleTlenValue) isize.increment(tLen);
			
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
		return totalRead + duplicate.get() + unmapped.get() + nonCanonical.get() ;		
	}
		
	/**
	 * check all globle value and assign the sumamry value
	 * eg.  
	 * 	private long trimedBase = 0; 	
	 */
	public void preSummary(long max , long trimB, long trimR, long baseDup, long baseUnmap, long baseNonCan) {
		//maxReadLength is continuelly updated during parseRecord
				
		//number of reads excluds discarded one.
		long no = duplicate.get() + unmapped.get() + nonCanonical.get() ;
		for (int i = 1 ; i <= this.maxReadLength ; i++)			 
			no += readLength.get(i);		
		this.noOfCountedReads = no;
		
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
		no = 0;
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
		badReadStats( rgElement, node_nonCanonicalPair, nonCanonical.get(),  getnonCanonicalBase());
				
		long lostBase = getDuplicateBase() + getUnmappedBase()+ getnonCanonicalBase();
		
		if( !readGroupId.equals(QprofilerXmlUtils.All_READGROUP) ){
			lostBase += lostBaseStats( rgElement, node_trim, parseTrim(readLength ), this.maxBases);				
		}else{
			lostBase += this.trimedBase;
			badReadStats( rgElement, node_trim, this.trimedRead, this.trimedBase   );
		}
					
		lostBase += lostBaseStats( rgElement, node_softClip, softClip, this.maxBases );
		lostBase += lostBaseStats( rgElement, node_hardClip, hardClip, this.maxBases );
		lostBase += lostBaseStats( rgElement, node_overlap, overlapBase, this.maxBases );
		XmlUtils.addCommentChild(( Element )rgElement.getLastChild(), "Only count overlapped bases on strand for pairs which have a positive TLEN value.");
		
		
		//create node for overall
		Element overallEle = XmlUtils.createGroupNode(rgElement, QprofilerXmlUtils.overall );
		XmlUtils.outputValueNode( overallEle, "readMaxLength", this.maxReadLength  );
		
		//readCount
		String comment = sreadCount + ": includes duplicateReads, nonCanonicalPairs and unmappedReads but excludes discardedReads (failed, secondary and supplementary).";
		overallEle.appendChild( overallEle.getOwnerDocument().createComment(comment) );				
		XmlUtils.outputValueNode( overallEle, sreadCount,  this.noOfCountedReads );
		
		
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
		Element ele =  XmlUtils.createMetricsNode(parent, "pairs", pairNum.get());

		
		Element ele1 = XmlUtils.createGroupNode(ele, "unPaired");
		XmlUtils.outputValueNode( ele1, "mateUnmappedPair", mateUnmapped.get() );
		XmlUtils.outputValueNode( ele1, "mateDifferentReferencePair", diffRef.get() );
		
		f5f3.toXml( ele );
		f3f5.toXml( ele );
		inward.toXml( ele );
		outward.toXml( ele );		

		//output tlen refer to method lostBaseStats(...)
		long sum = 0, no = 0;
		for(int i = 1; i < middleTlenValue; i ++ ){
			sum += isize.get(i) * i;
			no += isize.get(i);
		}		
		int mean = (no == 0) ? 0: (int) (sum / no);		
		int min = 0; //find the smallest non-zero value;
		for(int i = 1; i < middleTlenValue; i ++) {
			if(isize.get(i) > 0){ 
				min  = i; break; 
			}
		}		 
		int medium = 0;  	
		// to avoid isize.get(0) >= 0 since 1(counts)/2== 0(counts/2) == 0
		for (int i = 1,  total1 = 0 ; i < middleTlenValue; i++) {
			if(( total1 += isize.get(i)) > no/2 ){  
				medium = i;   break; 
			}
		}
		int mode = 0; //mode is the number of read which length is most popular
		long highest = 0;
		double sd = 0;
		for (int i = 1 ; i < middleTlenValue ; i++) {	
			//the number of isize.get(i) pairs have same isize value 
			sd += Math.pow(( i - mean), 2) * isize.get(i) / no;
			if(isize.get(i) > highest){ highest = isize.get(i); mode = i;  } 
		}
		double standardDeviation = Math.sqrt(sd);
		
		//tlen section 
		ele1 = XmlUtils.createGroupNode(ele, "tlen");
		XmlUtils.outputValueNode(ele1, smax, this.max_isize.get());
		XmlUtils.outputValueNode(ele1, smin, min);		
		
		ele1.appendChild(ele1.getOwnerDocument().createComment("below stats is based on tlen value < " + middleTlenValue));
		XmlUtils.outputValueNode(ele1, smean, mean);
		XmlUtils.outputValueNode(ele1, smode, mode);	
		XmlUtils.outputValueNode(ele1, smedian, medium);
		XmlUtils.outputValueNode(ele1, pairCount, no);		
		XmlUtils.outputValueNode( ele1, stdDev, (int)standardDeviation);			
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