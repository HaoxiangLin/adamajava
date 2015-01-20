package au.edu.qimr.qannotate.modes;

import static org.junit.Assert.assertTrue;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;
import org.qcmg.common.vcf.VcfInfoFieldRecord;
import org.qcmg.common.vcf.VcfRecord;
import org.qcmg.common.vcf.header.VcfHeaderRecord;
import org.qcmg.common.vcf.header.VcfHeaderRecord.MetaType;
import org.qcmg.common.vcf.header.VcfHeaderUtils;
import org.qcmg.vcf.VCFFileReader;
import org.qcmg.vcf.VCFFileWriter;

import au.edu.qimr.qannotate.utils.SnpEffMafRecord;

public class Vcf2mafTest {
	 @BeforeClass
	public static void createInput() throws IOException{
		createVcf();
		
	}
 
	 
	 @Test 
	 public void converterTest() throws Exception{
		 
		 
		 	final SnpEffMafRecord Dmaf = new SnpEffMafRecord();
			Dmaf.setDefaultValue();
			
			final Vcf2maf v2m = new Vcf2maf(2,1);	
			final String[] parms = {"chrY","22012840",".","C","A",".","SBIAS","MR=15;NNS=13;FS=GTGATATTCCC;EFF=sequence_feature[compositionally_biased_region:Glu/Lys-rich](LOW|||c.1252G>C|591|CCDC148|protein_coding|CODING|ENST00000283233|10|1),splice_acceptor_variant(HIGH|||n.356G>C||CCDC148-AS1|antisense|NON_CODING|ENST00000412781|5|1)","GT:GD:AC","0/0:C/A:A1[5],0[0],C6[6.67],0[0],T1[6],21[32.81]","0/1:C/A:C8[7.62],2[2],A2[8],28[31.18]"};
	 		final VcfRecord vcf = new VcfRecord(parms);
	 		final SnpEffMafRecord maf = v2m.converter(vcf);
	 		final String eff = new VcfInfoFieldRecord(vcf.getInfo()).getfield(VcfHeaderUtils.INFO_EFFECT);
	 			
	 		//select the annotation with "HIGH" impact
	 		//str: HIGH|||n.356G>C||CCDC148-AS1|antisense|NON_CODING|ENST00000412781|5|1
	 		//array: 0| 1|2|3     |4|5         |6        |7         |8              |9|10	 
	 		assertTrue(maf.getColumnValue(9).equals("HIGH" ));
	 		assertTrue(maf.getColumnValue(50).equals(Dmaf.getColumnValue(50) ));
	 		assertTrue(maf.getColumnValue(51).equals("n.356G>C" ));
	 		assertTrue(maf.getColumnValue(1).equals("CCDC148-AS1" ));
	 		assertTrue(maf.getColumnValue(53).equals("antisense" ));
	 		assertTrue(maf.getColumnValue(54).equals("NON_CODING"));
	 		assertTrue(maf.getColumnValue(49).equals("ENST00000412781" ));
	 		assertTrue(maf.getColumnValue(55).equals("5" ));
	 		assertTrue(maf.getColumnValue(56).equals("1" ));		
	 		
	 		//for other columns after A.M confirmation
	 		assertTrue(maf.getColumnValue(2).equals(Dmaf.getColumnValue(2) ));		
	 		assertTrue(maf.getColumnValue(3).equals(SnpEffMafRecord.Unknown ));		
	 		assertTrue(maf.getColumnValue(4).equals(Dmaf.getColumnValue(4) ));		
	 		assertTrue(maf.getColumnValue(5).equals(parms[0]));		
	 		assertTrue(maf.getColumnValue(6).equals(parms[1] ));		
	 		assertTrue(maf.getColumnValue(7).equals(parms[1] ));		
	 		assertTrue(maf.getColumnValue(8).equals(Dmaf.getColumnValue(8) ));		
	 		assertTrue(maf.getColumnValue(10).equals(Dmaf.getColumnValue(10) ));	
	 		assertTrue(maf.getColumnValue(11).equals(parms[3] ));		
	 		

	 		//check sample format column Vcf2maf(int tumour_column, int normal_column) = Vcf2maf(2,1); 
	 		assertTrue(maf.getColumnValue(36).equals(parms[9]  ));		//ND
	 		assertTrue(maf.getColumnValue(37).equals(parms[10] ));		//TD
	 		
	 		//check info field	 		
/*	 		maf[35] = Null; //ND
			maf[36] = Null; //TD
			assertTrue(maf.getColumnValue(12).equals("0" ));			
	 		assertTrue(maf.getColumnValue(13).equals("0" ));		
	 		assertTrue(maf.getColumnValue(14).equals(parms[2] ));	//dbsnp		
	 		
	 		assertTrue(maf.getColumnValue(11).equals("0" ));		
	 		assertTrue(maf.getColumnValue(12).equals("0" ));			
	 		assertTrue(maf.getColumnValue(11).equals("0" ));		
	 		assertTrue(maf.getColumnValue(12).equals("0" ));			
	 		assertTrue(maf.getColumnValue(11).equals("0" ));		
	 		assertTrue(maf.getColumnValue(12).equals("0" ));			
	 		assertTrue(maf.getColumnValue(11).equals("0" ));		
	 		assertTrue(maf.getColumnValue(12).equals("0" ));			
	 		assertTrue(maf.getColumnValue(11).equals("0" ));		
	 		assertTrue(maf.getColumnValue(12).equals("0" ));			
	 		assertTrue(maf.getColumnValue(11).equals("0" ));		
	 		assertTrue(maf.getColumnValue(12).equals("0" ));			
	 		
	 		assertTrue(maf.getColumnValue(9).equals("0" ));		
	 	*/	
	 		
	 }

	
	 @Test
	 public  void createMaf() throws IOException, Exception{
 		
         	
	 }
	 @Test
	 public  void singleSampleTest() throws IOException, Exception{
 		 
		 final File tmp = new File( DbsnpModeTest.inputName + ".tmp");
		 //create vcf by removing last column
		 try(final VCFFileReader reader = new VCFFileReader(new File( DbsnpModeTest.inputName)); 
				 final VCFFileWriter writer = new VCFFileWriter( tmp )){
				 for(final VcfHeaderRecord record: reader.getHeader()){ 
					 final String str = record.toString();
					 if(record.getMetaType().equals(MetaType.CHROM))
						 writer.addHeader(str.substring(0, str.lastIndexOf("\t") ));
					 else
						 writer.addHeader(str);
				 }
				 
				 for (final VcfRecord vcf : reader) {
					 final List<String> fields = vcf.getFormatFields();
					 fields.remove(2);
					 vcf.setFormatFields(fields);  
				//	 System.out.println(vcf.getFormatFields().size() + "\n" + vcf.getFormatFieldStrings()  ); 
					 writer.add(vcf);
				 }
	        } 
		 
		 
		final Vcf2maf v2m = new Vcf2maf(1,1);		 
		try(VCFFileReader reader = new VCFFileReader(tmp); ){
	//		System.out.println("SnpEff: " + SnpEffMafRecord.getSnpEffMafHeaderline());
			final SnpEffMafRecord[] maf = new SnpEffMafRecord[3];
			int i=0;
        	for (final VcfRecord vcf : reader){  		
        		maf[i++] = v2m.converter(vcf);
        		assertTrue( maf[i-1].getColumnValue(36).equals(maf[i-1].getColumnValue(37)) );
        	}	
        }
       
	 }	 
	 
	 //test get Sample column number
	
	
	public static void createVcf() throws IOException{
        final List<String> data = new ArrayList<String>();
        data.add("##fileformat=VCFv4.0");
        data.add(VcfHeaderUtils.STANDARD_FINAL_HEADER_LINE + "\tFORMAT\tCONTROL\tTEST");
       // data.add("chrY\t22012840\t.\tC\tA\t.\tSBIAS\tMR=15;NNS=13;FS=GTGATATTCCC;EFF=sequence_feature[compositionally_biased_region:Glu/Lys-rich](LOW|||c.1252G>C|591|CCDC148|protein_coding|CODING|ENST00000283233|10|1),splice_acceptor_variant(HIGH|||n.356G>C||CCDC148-AS1|antisense|NON_CODING|ENST00000412781|5|1)\tGT:GD:AC\t0/0:C/A:A1[5],0[0],C6[6.67],0[0],T1[6],21[32.81]\t0/1:C/A:C8[7.62],2[2],A2[8],28[31.18]");        
       // data.add("chr2\t91888700\trs78405093\tG\tA\t.\tPASS\tMR=1217;NNS=1;FS=TGAGCACCTAC;GMAF=0.113802559414991;EFF=intron_variant(MODIFIER|||n.376-566C>T||AC027612.3|processed_transcript|NON_CODING|ENST00000436174|4|1),intron_variant(MODIFIER|||n.478-123C>T||AC027612.3|transcribed_unprocessed_pseudogene|NON_CODING|ENST00000445955|4|1)\tGT:GD:AC\t1/1:G/T:G0[0],1217[35.41],T0[0],5[26.8],A0[0],7[34],C0[0],7[31]\t1/1:A/A:A0[0],1217[35.41],C0[0],5[26.8],G0[0],7[34],T0[0],7[31]"); 
    //    data.add("chrY\t2675825\t.\tTTG\tTCA\t.\tMIN;MIUN\tSOMATIC;END=2675826\tACCS\tTTG,5,37,TCA,0,2\tTAA,1,1,TCA,4,1,TCT,3,1,TTA,11,76,TTG,2,2,_CA,0,3,TTG,0,1");
    
        data.add("GL000236.1\t7127\t.\tT\tC\t.\tMR;MIUN\tSOMATIC;MR=4;NNS=4;FS=CCAGCCTATTT;EFF=non_coding_exon_variant(MODIFIER|||n.1313T>C||CU179654.1|processed_pseudogene|NON_CODING|ENST00000400789|1|1);CONF=ZERO\tGT:GD:AC\t0/0:T/T:T9[37.11],18[38.33]\t0/1:C/T:C1[12],3[41],T19[35.58],30[33.63]");
           try(BufferedWriter out = new BufferedWriter(new FileWriter(DbsnpModeTest.inputName));) {          
            for (final String line : data)   out.write(line + "\n");                  
         }  
	}

}
