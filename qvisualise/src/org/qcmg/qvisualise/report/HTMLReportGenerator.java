/**
 * © Copyright The University of Queensland 2010-2014.
 * © Copyright QIMR Berghofer Medical Research Institute 2014-2016.
 *
 * This code is released under the terms outlined in the included LICENSE file.
 */
package org.qcmg.qvisualise.report;

import java.util.ArrayList;
import java.util.List;

import org.qcmg.common.string.StringUtils;
import org.qcmg.qvisualise.ChartTab;

public class HTMLReportGenerator {
	
	private static final String END_DIV = "\n</div>";
	
	private List<Report> reports = new ArrayList<Report>();
	
	public HTMLReportGenerator(List<Report> reports) {
		this.reports = reports;
	}
	
	private String getReportHeaderInfo(Report report) {
		final StringBuilder sb = new StringBuilder();
		
		for (ChartTab tab : report.getTabs()) {
			if (tab.getChildren().isEmpty() && null != tab.getData() &&  null != tab.getChartInfo()) {
				sb.append(tab.getData());
				sb.append(tab.getChartInfo());
			} else {
				for (ChartTab child : tab.getChildren()) {
					if (null != child.getData() &&  null != child.getChartInfo() ) {
						sb.append(child.getData());
						sb.append(child.getChartInfo());
					}
				}
			}
		}
		
		return sb.toString();
	}
	
	public String generate() {
		StringBuilder sb = new StringBuilder();
		
		// get the header info
		sb.append("<html>\n<head>\n");
		sb.append(HTMLReportUtils.SS);
		sb.append(HTMLReportUtils.getStyle());
		sb.append("<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/1.4.4/jquery.min.js\"></script>\n");
		sb.append("<script type=\"text/javascript\"> $(document).ready(function () {\n");
		
		// loop through all tabs, if we have descriptions set, disable them initially so they are not 'on' by default
		for (Report report : reports) {
			for (ChartTab tab: report.getTabs()) {
				//special case for summary table description, it is not belong to tab
				if (null != tab.getDescription() || (tab.getName() != null && tab.getName().equals("summ1"))) {
					sb.append(" $(\"#").append(tab.getName()).append("1Desc_div\").toggle(false);\n");
					sb.append(" $(\"#").append(tab.getName()).append("2Desc_div\").toggle(false);\n");
				}
				for (ChartTab child: tab.getChildren()) {				
					if (null != child.getDescription() ) {
						sb.append(" $(\"#").append(child.getName()).append("Desc_div\").toggle(false);\n");
					} 
				}
			}
		}
		
		sb.append("});\n function toggleDiv(divId) { $(\"#\"+divId).toggle('fast'); }");
		sb.append("</script></head>\n");		
		sb.append("<body>");
			
		if (reports.size() > 1) {
			// don't want tab for report if it is the only one
			sb.append("\n<ul class=\"tabs\">");
			for (Report report : reports) {
				sb.append("<li><a href=\"#\">")
				.append(report.getType().getReportName())
				.append("</a></li>");
			}
			sb.append("\n</ul>");
		}
		
		for (Report report : reports) {
			sb.append(getReportBodyInfo(report));
		}
		// add google chart info
		HTMLReportUtils.generateHTMLHeader(sb);
		
		// get header info from each report
		for (Report report : reports) {
			sb.append(getReportHeaderInfo(report));
		}
		// end of chart info
		sb.append("}");
		// end of scripts
		sb.append( "\n</script>\n");
		sb.append("<!-- This JavaScript snippet activates those tabs -->" + 
				"\n<script> " +
				
				// perform JavaScript after the document is scriptable.
				"$(function() {" +
				"	$(\"ul.tabs\").tabs(\"> .pane\");" +
				"});" +
		"</script>");
		
		// end of file!!
		sb.append("\n</body>\n</html>");
		
		return sb.toString();
	}


	private String  getReportBodyInfo(Report report) {
		final StringBuilder sb = new StringBuilder();
		
		if (reports.size() > 1)
			sb.append("\n<div class=\"pane\"> ");
		
		
		sb.append("\n<div class=\"header\">File:").append(report.getFileName()).append(END_DIV);		
		if( ! StringUtils.isNullOrEmpty( report.getRecordParsed() ))
			sb.append(String.format("\n<div class=\"header\">RecordsParsed:%s; RunBy:%s; RunOn:%s; Version:qprofiler-%s",report.getRecordParsed(), report.getRunBy(), report.getRunOn(), report.getVersion())).append(END_DIV);		
		else
			sb.append(String.format("\n<div class=\"header\">RunBy:%s; RunOn:%s; Version:qprofiler-%s", report.getRunBy(), report.getRunOn(), report.getVersion())).append(END_DIV);		
			
		
		sb.append(END_DIV);	
		sb.append("\n<div class=\"header\">&nbsp</div>");		
		// and add the tabs to the body
		sb.append("\n<ul class=\"tabs\">");
		// display parent tabs
		for (ChartTab tab : report.getTabs()) {  
			sb.append("\n<li><a href=\"#\">");
			sb.append(tab.getTitle()).append("</a></li>");
		}
		sb.append("\n</ul>");
		
		
		for (ChartTab tab : report.getTabs()) {			
			if (null == tab.getChildren() || tab.getChildren().isEmpty()) {
				if (null == tab.getRenderingInfo()) {
					if (null == tab.getDescription()) {
 						sb.append("\n<div class=\"pane\" id=\"").append(tab.getName()).append("Chart_div\">");
 						sb.append(END_DIV);						
					} else {
						sb.append("\n<div class=\"pane\">");
						sb.append("<p id=\"").append(tab.getName()).append( "Chart_div\"></p>");
						sb.append(END_DIV);
					}
				} else {
					sb.append(tab.getRenderingInfo());
				}
			} else {
				sb.append("\n<div class=\"pane\">"); 
				sb.append("\n<ul class=\"tabs\">");
				for (ChartTab child : tab.getChildren()) {
					sb.append("\n<li><a href=\"#\">");
					sb.append(child.getTitle()).append("</a></li>");
				}
				sb.append("\n</ul>");
				for (ChartTab child : tab.getChildren()) {
					if (null != child.getRenderingInfo()) {
						sb.append(child.getRenderingInfo());  //xu code:: now child have rendering 
					}else if (null == child.getDescription()) {
						sb.append("\n<div class=\"pane\" id=\"").append(child.getName()).append("Chart_div\">");
						sb.append(END_DIV);
					} else {
						sb.append("\n<div class=\"pane\">");
						sb.append(HTMLReportUtils.generateDescriptionButton(child.getName(), child.getDescription(), null));							
						sb.append("<p id=\"").append(child.getName()).append("Chart_div\"></p>");
						sb.append(END_DIV);
					}
				}
				sb.append(END_DIV);
			}
		}
		
		if (reports.size() > 1)
			sb.append(END_DIV).append("\n");
		
		return sb.toString();
	}

}
