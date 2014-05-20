package org.qcmg.snp;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import junit.framework.Assert;

import org.junit.Test;
import org.junit.rules.ExpectedException;
import org.junit.rules.TemporaryFolder;
import org.qcmg.common.commandline.Executor;
import org.qcmg.common.model.PileupElement;
import org.qcmg.common.model.Rule;
import org.qcmg.common.model.VCFRecord;
import org.qcmg.tab.TabbedFileReader;
import org.qcmg.tab.TabbedRecord;
import org.qcmg.vcf.VCFFileReader;

public class VcfPipelineTest {
	
	private final static char NL = '\n';
	
	@org.junit.Rule
	public  TemporaryFolder testFolder = new TemporaryFolder();
	@org.junit.Rule
    public ExpectedException thrown= ExpectedException.none();
	
	
	
	
	
	private void createVcfFile(File vcfFile, List<String> data) throws IOException {
		FileWriter writer = new FileWriter(vcfFile);
		try {
			// add data
			for (String s : data) {
				writer.write(s + NL);
			}
		} finally {
			writer.close();
		}
	}
	
	@Test
	public void testPileupPipelineEmptyPileupFile() throws Exception {
		File logFile = testFolder.newFile("qsnp.log");
		File iniFile = testFolder.newFile("qsnp.ini");
		IniFileGenerator.createRulesOnlyIni(iniFile);
		
		File pileupInput = testFolder.newFile("input.pileup");
		File vcfOutput = testFolder.newFile("output.vcf");
		
//		PileupFileGenerator.createPileupFile(pileupInput);
		
		IniFileGenerator.addInputFiles(iniFile, false, "pileup = " + pileupInput.getAbsolutePath());
		IniFileGenerator.addOutputFiles(iniFile, false, "vcf = " + vcfOutput.getAbsolutePath());
		
		// that should be it
		String command = "-log " + logFile.getAbsolutePath() + " -i " + iniFile.getAbsolutePath();
		Executor exec = new Executor(command, "org.qcmg.snp.Main");
		assertEquals(1, exec.getErrCode());
	}
	
	@Test
	public void testPileupPipelineGenerateVCFOnly() throws Exception {
		File logFile = testFolder.newFile("qsnp.log");
		File iniFile = testFolder.newFile("qsnp.ini");
		IniFileGenerator.createRulesOnlyIni(iniFile);
		
		File pileupInput = testFolder.newFile("input.pileup");
		File vcfOutput = testFolder.newFile("output.vcf");
		
		PileupFileGenerator.createPileupFile(pileupInput);
		
		IniFileGenerator.addInputFiles(iniFile, false, "pileup = " + pileupInput.getAbsolutePath());
		IniFileGenerator.addOutputFiles(iniFile, false, "vcf = " + vcfOutput.getAbsolutePath());
		
		// add runType to ini file
		IniFileGenerator.addStringToIniFile(iniFile, "[parameters]\nrunMode = pileup", true);	// append to file
		
		// that should be it
		ExpectedException.none();
		String command = "-log " + logFile.getAbsolutePath() + " -i " + iniFile.getAbsolutePath();
		Executor exec = new Executor(command, "org.qcmg.snp.Main");
		assertEquals(0, exec.getErrCode());
		assertTrue(0 == exec.getOutputStreamConsumer().getLines().length);
		
		// check the vcf output file
		assertEquals(1, noOfLinesInVCFOutputFile(vcfOutput));
	}
	
	@Test
	public void testPileupPipelineGenerateVCFOnlyIncludeIndels() throws Exception {
		File logFile = testFolder.newFile("qsnp.log");
		File iniFile = testFolder.newFile("qsnp.ini");
		IniFileGenerator.createRulesOnlyIni(iniFile);
		
		File pileupInput = testFolder.newFile("input.pileup");
		File vcfOutput = testFolder.newFile("output.vcf");
		
		PileupFileGenerator.createPileupFile(pileupInput);
		
		IniFileGenerator.addInputFiles(iniFile, false, "pileup = " + pileupInput.getAbsolutePath());
		IniFileGenerator.addOutputFiles(iniFile, false, "vcf = " + vcfOutput.getAbsolutePath());
		IniFileGenerator.addStringToIniFile(iniFile, "[parameters]\nincludeIndels = true", true);
		
		// add runType to ini file
		IniFileGenerator.addStringToIniFile(iniFile, "\nrunMode = pileup", true);	// append to file
		
		// that should be it
		ExpectedException.none();
		String command = "-log " + logFile.getAbsolutePath() + " -i " + iniFile.getAbsolutePath();
		Executor exec = new Executor(command, "org.qcmg.snp.Main");
		assertEquals(0, exec.getErrCode());
		assertTrue(0 == exec.getOutputStreamConsumer().getLines().length);
		
		// check the vcf output file
		assertEquals(2, noOfLinesInVCFOutputFile(vcfOutput));
	}
	
	@Test
	public void testPileupPipelineDCCMode() throws Exception{
		File logFile = testFolder.newFile("qsnp.log");
		File iniFile = testFolder.newFile("qsnp.ini");
		IniFileGenerator.createRulesOnlyIni(iniFile);
		
		File pileupInput = testFolder.newFile("input.pileup");
		File vcfOutput = testFolder.newFile("output.vcf");
		File pileupOutput = testFolder.newFile("output.pileup");
		File dccSomaticOutput = testFolder.newFile("output.dcc.somatic");
		File dccGermlineOutput = testFolder.newFile("output.dcc.germline");
		File illuminaFileNormalAndTumour = testFolder.newFile("illumina.normal.tumour");
		File chrConv = testFolder.newFile("chr.conv");
		File germlineDB = testFolder.newFile("germline.DB");
		
		PileupFileGenerator.createPileupFile(pileupInput);
		IlluminaFileGenerator.createIlluminaFile(illuminaFileNormalAndTumour);
		
		IniFileGenerator.addInputFiles(iniFile, false, "pileup = " + pileupInput.getAbsolutePath()
				+ "\nilluminaNormal = " + illuminaFileNormalAndTumour.getAbsolutePath()
				+ "\nilluminaTumour = " + illuminaFileNormalAndTumour.getAbsolutePath()
				+ "\ngermlineDB = " + germlineDB.getAbsolutePath()
				+ "\nchrConv = " + chrConv.getAbsolutePath());
		
		IniFileGenerator.addOutputFiles(iniFile, false, "vcf = " + vcfOutput.getAbsolutePath() 
				+ "\ndccSomatic = " + dccSomaticOutput.getAbsolutePath()
				+ "\npileup = " + pileupOutput.getAbsolutePath()
				+ "\ndccGermline = " + dccGermlineOutput.getAbsolutePath());
		IniFileGenerator.addStringToIniFile(iniFile, "[parameters]\nincludeIndels = true", true);
		
		// add the annotate mode=dcc to the ini file
		IniFileGenerator.addStringToIniFile(iniFile, "\nannotateMode = dcc", true);
		// add runType to ini file
		IniFileGenerator.addStringToIniFile(iniFile, "\nrunMode = pileup", true);	// append to file
		
		// that should be it
		ExpectedException.none();
		String command = "-log " + logFile.getAbsolutePath() + " -i " + iniFile.getAbsolutePath();
		Executor exec = new Executor(command, "org.qcmg.snp.Main");
		assertEquals(0, exec.getErrCode());
		assertTrue(0 == exec.getOutputStreamConsumer().getLines().length);
		
		// check the vcf output file
		assertEquals(2, noOfLinesInVCFOutputFile(vcfOutput));
		// check the dcc somatic output file
		assertEquals(1, noOfLinesInDCCOutputFile(dccSomaticOutput));
		// check the dcc germline output file
		assertEquals(1, noOfLinesInDCCOutputFile(dccGermlineOutput));
	}
	
	private int noOfLinesInVCFOutputFile(File vcfOutput) throws Exception {
		VCFFileReader reader = new VCFFileReader(vcfOutput);
		int noOfLines = 0;
		try {
			for (VCFRecord vcf : reader) noOfLines++;
		} finally {
			reader.close();
		}
		return noOfLines;
	}
	
	private int noOfLinesInDCCOutputFile(File dccFile) throws Exception {
		TabbedFileReader reader = new TabbedFileReader(dccFile);
		int noOfLines = 0;
		try {
			for (TabbedRecord vcf : reader) {
				if (vcf.getData().startsWith("analysis")) continue;	// header line
				noOfLines++;
			}
		} finally {
			reader.close();
		}
		return noOfLines;
	}
	
	@Test
	public void testIsRecordAKeeper() {
		// arguments are.....
		// int variantCount, int coverage, Rule rule, List<PileupElement> baseCounts, double percentage
		
		Rule r = new Rule(0, 20, 1);
		List<PileupElement> pes = new ArrayList<PileupElement>();
		
		assertEquals(false, PileupPipeline.isPileupRecordAKeeper(0, 0, r, null, 0));
		assertEquals(false, PileupPipeline.isPileupRecordAKeeper(0, 0, r, pes, 0));
		
		try {	// variant count is greater than total count
			assertEquals(false, PileupPipeline.isPileupRecordAKeeper(1, 0, r, pes, 0));
			Assert.fail("Should have thrown an IllegalArgumentException");
		} catch (IllegalArgumentException e) {}
		
		assertEquals(false, PileupPipeline.isPileupRecordAKeeper(1, 1, r, pes, 0));
		PileupElement pe = new PileupElement('A');
		pe.incrementForwardCount((byte)'I');
		pes.add(pe);
		assertEquals(true, PileupPipeline.isPileupRecordAKeeper(1, 1, r, pes, 0));
		
		//. change rule
		assertEquals(false, PileupPipeline.isPileupRecordAKeeper(1, 1, new Rule(0,20,2), pes, 0));
		pe.incrementReverseCount((byte)'I');
		assertEquals(true, PileupPipeline.isPileupRecordAKeeper(2, 2, new Rule(0,20,2), pes, 0));
		assertEquals(false, PileupPipeline.isPileupRecordAKeeper(2, 2, new Rule(0,20,4), pes, 0));
		
		// only use percentage if we are dealing with the upper bounded rule
		assertEquals(true, PileupPipeline.isPileupRecordAKeeper(2, 100, new Rule(0,Integer.MAX_VALUE,4), pes, 0));
		assertEquals(false, PileupPipeline.isPileupRecordAKeeper(1, 100, new Rule(0,Integer.MAX_VALUE,4), pes, 0));
		
		PileupElement pe2 = new PileupElement('.');
		pe2.incrementForwardCount((byte)'I');
		pe2.incrementForwardCount((byte)'I');
		pe2.incrementForwardCount((byte)'I');
		pe2.incrementForwardCount((byte)'I');
		pe2.incrementForwardCount((byte)'I');
		pe2.incrementForwardCount((byte)'I');
		pe2.incrementForwardCount((byte)'I');
		pes.add(pe2);
		assertEquals(false, PileupPipeline.isPileupRecordAKeeper(2, 9, new Rule(0,20,2), pes, 50));
		assertEquals(true, PileupPipeline.isPileupRecordAKeeper(2, 9, new Rule(0,20,2), pes, 10));
	}
	
	@Test
	public void testIsVariantOnBothStrands() {
		assertEquals(false, PileupPipeline.isVariantOnBothStrands(null));
		List<PileupElement> pes = new ArrayList<PileupElement>();
		assertEquals(false, PileupPipeline.isVariantOnBothStrands(pes));
		PileupElement pe = new PileupElement('.');
		pe.incrementForwardCount((byte)'I');
		pes.add(pe);
		assertEquals(false, PileupPipeline.isVariantOnBothStrands(pes));
		PileupElement pe2 = new PileupElement('C');
		pe2.incrementForwardCount((byte)'I');
		pes.add(pe2);
		assertEquals(false, PileupPipeline.isVariantOnBothStrands(pes));
		pe2.incrementReverseCount((byte)'?');
		assertEquals(true, PileupPipeline.isVariantOnBothStrands(pes));
		
		pes.clear();
		PileupElement pe3 = new PileupElement('T');
		pe3.incrementReverseCount((byte)'I');
		pes.add(pe3);
		assertEquals(false, PileupPipeline.isVariantOnBothStrands(pes));
		pe3.incrementForwardCount((byte)'I');
		assertEquals(true, PileupPipeline.isVariantOnBothStrands(pes));
	}
}