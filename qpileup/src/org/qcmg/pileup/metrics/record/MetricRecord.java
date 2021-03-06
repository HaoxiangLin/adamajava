/**
 * © Copyright The University of Queensland 2010-2014.  This code is released under the terms outlined in the included LICENSE file.
 */
package org.qcmg.pileup.metrics.record;


import java.text.DecimalFormat;

import org.qcmg.pileup.PileupConstants;
import org.qcmg.pileup.PileupUtil;

public class MetricRecord {
	
	protected final Integer position;
	protected Long count;
	protected final String type;
	protected final String chromosome;
	protected int totalReads;
	protected int endPosition;

	public MetricRecord(String type, String chromosome, Integer position, long count, int totalReads) {
		this.type = type;
		this.chromosome = chromosome;
		this.position = position;
		this.endPosition = position;
		this.count = count;
		this.totalReads = totalReads;
	}
	
	public int getTotalReads() {
		return totalReads;
	}

	public Integer getPosition() {
		return position;
	}

	public Long getCount() {
		return count;
	}

	public void setCount(Long count) {
		this.count = count;
	}

	public void setTotalReads(int totalReads) {
		this.totalReads = totalReads;		
	}	

	public int getEndPosition() {
		return endPosition;
	}

	public void setEndPosition(Integer endPosition) {
		this.endPosition = endPosition;
	}

	private String getGFFType() {
		if (type.equals(PileupConstants.METRIC_CLIP)) {
			return "CLIP";
		}
		if (type.equals(PileupConstants.METRIC_NONREFBASE)) {
			return "NR";
		}
		return "";
	}

	public String getType() {
		return this.type;
	}

	public double getRegularityScore() {
		if (PileupUtil.isRegularityType(type)) {
			double percent = getPercentage(count.longValue(), totalReads);
			return percent * percent;
		} 
		return 0;
	}
	
	public static double getPercentage(long count, long total) {
		if (count > total) {
			return 100;
		} else if (count == 0 || total == 0) {
			return 0;
		} else {
			return ((double)count/total) * 100;
		}	
	}

	public String toTmpString() {
		return chromosome + "\t" + position + "\t" + endPosition + "\t" + type + "\t" + count + "\t" + totalReads + "\n";
	}

	public String toGFFString() {
		double perc = getPercentage(count.longValue(), totalReads);
		DecimalFormat f = new DecimalFormat("##.00");
		String result = chromosome + "\t";
		result += "qpileup" + "\t";
		result += "." + "\t";
		result += position + "\t";
		result += endPosition + "\t";
		result += f.format(perc) + "\t";
		result += "." + "\t";
		result += "." + "\t";
		result += "Name=" + getGFFType() + ";color=" + "#C0C0C0" + ";PercentScore=" + f.format(perc) + "\n";
		return result;
	}	

	@Override
	public String toString() {
		return this.type + "\t" + chromosome + " " + position + " " + count;
	}

	public boolean hasStrandBias() {
		return false;
	}
}
