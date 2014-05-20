/**
 * © Copyright The University of Queensland 2010-2014.  This code is released under the terms outlined in the included LICENSE file.
 */
package org.qcmg.qsv.discordantpair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.qcmg.qsv.util.QSVConstants;



public class QPrimerCategory {

	private Map<String, Integer> map;
	private String primaryCategoryNo;
	private String mixedCategories;
	private int startLeft;
	private int endLeft;
	private int startRight;
	private int endRight;
	private String reverseFlag;
	private String zpType;
	private String leftChr;
	private String rightChr;
	private String id;
	private String outLeft;
	private String outRight;
	private String pairType;

	public QPrimerCategory(String zpType, String leftChr, String rightChr, String id, String pairType) {
		this.map = new HashMap<String,Integer>();
		this.zpType = zpType;
		this.leftChr = leftChr;
		this.rightChr = rightChr;
		this.mixedCategories = "";
		this.id = id;
		this.outLeft = new String();
		this.outRight = new String();
		this.pairType = pairType;
	}

	public Map<String, Integer> getMap() {
		return map;
	}

	public void setMap(Map<String, Integer> map) {
		this.map = map;
	}


	public String getPrimaryCategoryNo() {
		return primaryCategoryNo;
	}


	public void setPrimaryCategoryNo(String primaryCategoryNo) {
		this.primaryCategoryNo = primaryCategoryNo;
	}


	public String getMixedCategories() {
		return mixedCategories;
	}


	public void setMixedCategories(String mixedCategories) {
		this.mixedCategories = mixedCategories;
	}


	public int getStartLeft() {
		return startLeft;
	}


	public void setStartLeft(int startLeft) {
		this.startLeft = startLeft;
	}


	public int getEndLeft() {
		return endLeft;
	}


	public void setEndLeft(int endLeft) {
		this.endLeft = endLeft;
	}


	public int getStartRight() {
		return startRight;
	}


	public void setStartRight(int startRight) {
		this.startRight = startRight;
	}


	public int getEndRight() {
		return endRight;
	}


	public void setEndRight(int endRight) {
		this.endRight = endRight;
	}


	public String getReverseFlag() {
		return reverseFlag;
	}


	public void setReverseFlag(String reverseFlag) {
		this.reverseFlag = reverseFlag;
	}


	public String getZpType() {
		return zpType;
	}


	public void setZpType(String zpType) {
		this.zpType = zpType;
	}


	public String getLeftChr() {
		return leftChr;
	}


	public void setLeftChr(String leftChr) {
		this.leftChr = leftChr;
	}


	public String getRightChr() {
		return rightChr;
	}


	public void setRightChr(String rightChr) {
		this.rightChr = rightChr;
	}


	public String getId() {
		return id;
	}


	public void setId(String id) {
		this.id = id;
	}


	public String getOutLeft() {
		return outLeft;
	}


	public void setOutLeft(String outLeft) {
		this.outLeft = outLeft;
	}


	public String getOutRight() {
		return outRight;
	}


	public void setOutRight(String outRight) {
		this.outRight = outRight;
	}


	public void findClusterCategory(List<MatePair> clusterMatePairs, int clusterLeftStart, int clusterLeftEnd, int clusterRightStart, int clusterRightEnd) throws Exception {
		//find the categories of each read pair		
		countCategories(clusterMatePairs);
		
		if (map.size() == 0) {
			primaryCategoryNo= "unknown";
		} else {
 		
			findCategoryNo();
		
			findQPrimerSites(clusterLeftStart, clusterLeftEnd, clusterRightStart, clusterRightEnd);
		}
	}

	void countCategories(List<MatePair> clusterMatePairs) {
		String cat = null;
		if (pairType.equals("lmp")) {
			for (MatePair p: clusterMatePairs) {	
				cat = p.getSVCategoryForLMP();				
			}
		} else if (pairType.equals("pe")) {
			for (MatePair p: clusterMatePairs) {
				cat = p.getSVCategoryForPE();	
			}
		} else if (pairType.equals("imp")){
			for (MatePair p: clusterMatePairs) {
				cat = p.getSVCategoryForIMP();	
			}
		}
		
		if (cat != null) {
			addToCategoryMap(cat);
		}
	}
	

	public void findQPrimerSites(int clusterLeftStart, int clusterLeftEnd, int clusterRightStart, int clusterRightEnd) {
		if (primaryCategoryNo.equals(QSVConstants.ORIENTATION_1) || primaryCategoryNo.equals(QSVConstants.ORIENTATION_3)) {
			
			if (pairType.equals("lmp")) {
				setStandardEnds(clusterLeftStart, clusterLeftEnd, clusterRightStart, clusterRightEnd);		
			} else {
				if (primaryCategoryNo.equals(QSVConstants.ORIENTATION_1)) {
					setCat1PEEnds(clusterLeftStart, clusterLeftEnd, clusterRightStart, clusterRightEnd);
				} else {
					setCat3PEEnds(clusterLeftStart, clusterLeftEnd, clusterRightStart, clusterRightEnd);
				}
			}
			
			if (primaryCategoryNo.equals(QSVConstants.ORIENTATION_3)) {
				reverseFlag = "true";
			} else if (primaryCategoryNo.equals(QSVConstants.ORIENTATION_4)) {
				
			} else {
				reverseFlag = "false";
			}
		} else if (primaryCategoryNo.equals(QSVConstants.ORIENTATION_2)) {
			if (pairType.equals("lmp")) {
				setSwappedEnds(clusterLeftStart, clusterLeftEnd, clusterRightStart, clusterRightEnd);		
			} else {
				setCat2SwappedEnds(clusterLeftStart, clusterLeftEnd, clusterRightStart, clusterRightEnd);	
			}
				
			reverseFlag = "false";
		} else if (primaryCategoryNo.equals(QSVConstants.ORIENTATION_4)) {
			if (pairType.equals("lmp")) {
				setCat4Ends(clusterLeftStart, clusterLeftEnd, clusterRightStart, clusterRightEnd);	
			} else {
				setCat4PEEnds(clusterLeftStart, clusterLeftEnd, clusterRightStart, clusterRightEnd);	
			}
			
			reverseFlag = "lefttrue";
		} else if (primaryCategoryNo.equals(QSVConstants.ORIENTATION_5)) {
			setCat5SwappedEnds(clusterLeftStart, clusterLeftEnd, clusterRightStart, clusterRightEnd);	
			reverseFlag = "false";
		}
	}

	private void setCat1PEEnds(int clusterLeftStart, int clusterLeftEnd,
			int clusterRightStart, int clusterRightEnd) {
		startLeft = clusterLeftStart;
		endLeft = clusterLeftEnd - 50;
		startRight = clusterRightStart + 50;
		endRight = clusterRightEnd;		
	}
	
	private void setCat3PEEnds(int clusterLeftStart, int clusterLeftEnd, int clusterRightStart, int clusterRightEnd) {
		startLeft = clusterLeftStart;
		endLeft = clusterLeftEnd - 50;
		startRight = clusterRightStart;
		endRight = clusterRightEnd - 50;	
		String chr = leftChr;
		leftChr = rightChr;
		rightChr = chr;
	}
	
	private void setCat4PEEnds(int clusterLeftStart, int clusterLeftEnd, int clusterRightStart, int clusterRightEnd) {
		startLeft = clusterLeftStart+50;
		endLeft = clusterLeftEnd;
		startRight = clusterRightStart+50;
		endRight = clusterRightEnd;	
	}
	
	private void setCat2SwappedEnds(int clusterLeftStart, int clusterLeftEnd, int clusterRightStart, int clusterRightEnd) {
		startLeft = clusterRightStart;
		endLeft = clusterRightEnd - 50;		
		startRight = clusterLeftStart + 50;
		endRight = clusterLeftEnd;
		String chr = leftChr;
		leftChr = rightChr;
		rightChr = chr;
	}

	private void setCat5SwappedEnds(int clusterLeftStart, int clusterLeftEnd, int clusterRightStart, int clusterRightEnd) {
		int leftPrimerStart = clusterLeftStart + 50;
		int rightPrimerEnd = clusterRightEnd - 50;
		Double middle = Math.ceil(((new Double(clusterLeftStart) + new Double(clusterRightEnd))/ new Double(2)));
		
		int leftPrimerEnd = middle.intValue();
		int rightPrimerStart = middle.intValue();
		
		//swap!
		startLeft = rightPrimerStart;
		endLeft = rightPrimerEnd;
		startRight = leftPrimerStart;
		endRight = leftPrimerEnd;
		
		String chr = leftChr;
		leftChr = rightChr;
		rightChr = chr;
	}

	private void setStandardEnds(int clusterLeftStart, int clusterLeftEnd, int clusterRightStart, int clusterRightEnd) {
		
		Double endL = Math.ceil(((new Double(clusterLeftEnd) + new Double(clusterLeftStart))/ new Double(2)));		
		endLeft = endL.intValue();
		startLeft = endLeft-499;
		if (startLeft<clusterLeftStart) {
			startLeft = clusterLeftStart;
		}
		
		Double startR =  Math.ceil(((new Double(clusterRightEnd) + new Double(clusterRightStart))/2));		
		startRight = startR.intValue();
		
		endRight = startRight + 499;
		if (endRight>clusterRightEnd) {
			endRight = clusterRightEnd;
		}		
	}
	
	private void setCat4Ends(int clusterLeftStart, int clusterLeftEnd, int clusterRightStart, int clusterRightEnd) {
		
		Double startL = Math.ceil(((new Double(clusterLeftEnd) + new Double(clusterLeftStart))/2));		
		startLeft = startL.intValue();
		
		endLeft = startLeft+499;
		if (endLeft>clusterLeftEnd) {
			endLeft = clusterLeftEnd;
		}
		
		Double endR = Math.ceil(((new Double(clusterRightEnd) + new Double(clusterRightStart))/2));		
		endRight = endR.intValue();
		
		startRight = endRight - 499;
		if (startRight<clusterRightStart) {
			startRight = clusterRightStart;
		}		
	}
	
	private void setSwappedEnds(int clusterLeftStart, int clusterLeftEnd, int clusterRightStart, int clusterRightEnd) {
		//swap the left and right chr
		String tempRight = leftChr;
		String tempLeft = rightChr;
		
		leftChr = tempLeft;
		rightChr = tempRight;
		Double endL = Math.ceil(((new Double(clusterRightEnd) + new Double(clusterRightStart))/2));			
		endLeft = endL.intValue();
		
		startLeft = endLeft-499;
		if (startLeft<clusterRightStart) {
			startLeft = clusterRightStart;
		}
		
		Double startR = Math.ceil(((new Double(clusterLeftEnd) + new Double(clusterLeftStart))/2));
		startRight = startR.intValue();

		endRight = startRight + 499;
		if (endRight>clusterLeftEnd) {
			endRight = clusterLeftEnd;
		}
	}

	private void addToCategoryMap(String key) {
		if (map.containsKey(key)) {
		 	Integer newValue = new Integer(map.get(key).intValue() + 1);
		 	map.put(key, newValue);
		} else {
			map.put(key, new Integer(1));
		}		
	}
	
	public void findCategoryNo() throws Exception {
		this.primaryCategoryNo = "";
			if (map.size() == 0) {
				throw new Exception ("Pairs could not be assigned to qPrimer category");
			} else {
				if (map.size() == 1) {
					for (Entry<String, Integer> entry: map.entrySet()) {
						primaryCategoryNo = entry.getKey();
					}    			
				} else {
					List<String> categoryNos = new ArrayList<String>();
					int max = -1;
					//find the max no of pairs
					for (Entry<String, Integer> entry: map.entrySet()) {
						mixedCategories += "Cat" + entry.getKey() + "(" + entry.getValue() + "),";
						if (entry.getValue().intValue() > max) {
							max = entry.getValue().intValue();
						}
					}
					//see if it happens more than once
					for (Entry<String, Integer> entry: map.entrySet()) {
						if (entry.getValue().intValue() == max) {
							categoryNos.add(entry.getKey());
						}
					}			
					primaryCategoryNo = categoryNos.get(0);
				}
		}
	}
	
	public String toString(String svId) {
		String left = new String();
		String right = new String();
		
		left = leftChr + ":" + startLeft + "-" + endLeft;
		right = rightChr + ":" + startRight + "-" + endRight;

			return svId + "\t" + left + "\t" + right + "\t" + reverseFlag + "\t" + primaryCategoryNo + "\t" + mixedCategories;
	
	}
}