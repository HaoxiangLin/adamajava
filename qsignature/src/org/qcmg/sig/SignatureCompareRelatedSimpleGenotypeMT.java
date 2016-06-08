/**
 * ?? Copyright The University of Queensland 2010-2014.  This code is released under the terms outlined in the included LICENSE file.
 */
package org.qcmg.sig;

import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntShortHashMap;

import java.io.File;
import java.util.AbstractQueue;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.IntSummaryStatistics;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.qcmg.common.log.QLogger;
import org.qcmg.common.log.QLoggerFactory;
import org.qcmg.common.model.ChrPosition;
import org.qcmg.common.util.DonorUtils;
import org.qcmg.common.util.FileUtils;
import org.qcmg.common.util.LoadReferencedClasses;
import org.qcmg.sig.model.Comparison;
import org.qcmg.sig.util.ComparisonUtil;
import org.qcmg.sig.util.SignatureUtil;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * This class gets a list of all .qsig.vcf files from the supplied path.
 * It then performs a comparison between them all, regardless of whether they are bam or snp chip files
 * An xml output file is produced
 * If any comparison scores are less than the cutoff, they are added to a list, which is then emailed to interested parties informing them of the potential problem files
 *  
 * @author o.holmes
 *
 */
public class SignatureCompareRelatedSimpleGenotypeMT {
	
	private static QLogger logger;
	private int exitStatus;
	
	
	private float cutoff = 0.2f;
	private int minimumCoverage = 10;
	private final int nThreads = 4;
	
//	private final int cacheSize = 700;
	
	private String outputXml;
	private String [] paths;
	private String [] additionalSearchStrings;
	private String donor;
//	private static final String QSIG_SUFFIX = ".qsig.vcf";
	
	private String excludeVcfsFile;
	private List<String> excludes;
	private String logFile;
	
	private final Map<File, int[]> fileIdsAndCounts = new THashMap<>();
	private final List<Comparison> allComparisons = new ArrayList<>();
	
//	private final Map<File, Map<ChrPosition, float[]>> cache = new THashMap<>(cacheSize * 2);
	private final Map<File, TIntShortHashMap> cache = new THashMap<>();
	
	List<String> suspiciousResults = new ArrayList<String>();
	
	private int engage() throws Exception {
		
		// get excludes
		logger.info("Retrieving excludes list from: " + excludeVcfsFile);
		excludes = SignatureUtil.getEntriesFromExcludesFile(excludeVcfsFile);
		
		// get qsig vcf files for this donor
		logger.info("Retrieving qsig vcf files from: " + Arrays.stream(paths).collect(Collectors.joining(",")));
		Set<File> uniqueFiles = new HashSet<>();
		for (String path : paths) {
			uniqueFiles.addAll(FileUtils.findFilesEndingWithFilterNIO(path, SignatureUtil.QSIG_VCF));
		}
		List<File> files = new ArrayList<>(uniqueFiles);
		
		if (files.isEmpty()) {
			logger.warn("Didn't find any files ending with " + SignatureUtil.QSIG_VCF + " in " + Arrays.toString(paths));
			return 0;
		}
		
		logger.info("Total number of files to be compared: " + files.size());
		// remove excluded files
		files = SignatureUtil.removeExcludedFilesFromList(files, excludes);
		logger.info("Total number of files to be compared (minus excluded files): " + files.size());
		
		if (files.isEmpty()) {
			logger.warn("No files left after removing exlcuded files");
			return 0;
		}
		
		
		
		
		/*
		 * Match files on additionalSearchStrings
		 */
		if (null != additionalSearchStrings && additionalSearchStrings.length > 0) {
			Predicate<File> p = (File f) -> {
				return Arrays.stream(additionalSearchStrings).anyMatch(s -> f.getAbsolutePath().contains(s));
			};
			files = files.stream().filter(f -> p.test(f)).collect(Collectors.toList());
		}
		
		
		logger.info("Should have " + (files.size() -1) + " + " + (files.size() -2) + " ...  comparisons");
		
		Collections.sort(files, FileUtils.FILE_COMPARATOR);
		
		// add files to map
		addFilesToMap(files);
		
		if (donor == null) {
			donor = DonorUtils.getDonorFromFilename(files.get(0).getAbsolutePath());
			if (null == donor) {
				logger.warn("Could not get donor information from file: " + files.get(0).getAbsolutePath());
			}
		}
		
		StringBuilder donorSB = new StringBuilder(donor + "\n");
		
		
		
		populateCache(files);
		
		performComparisons(files);
		
		
//		int size = files.size();
//		
//		for (int i = 0 ; i < size -1 ; i++) {
//			
//			File f1 = files.get(i);
////			Map<ChrPosition, float[]> ratios1 = getSignatureData(f1);
//			TIntShortHashMap ratios1 = getSignatureData(f1);
//			
//			for (int j = i + 1 ; j < size ; j++) {
//				File f2 = files.get(j);
////				Map<ChrPosition, float[]> ratios2 = getSignatureData(f2);
//				TIntShortHashMap ratios2 = getSignatureData(f2);
//				
////				Comparison comp = QSigCompareDistance.compareRatiosFloat(ratios1, ratios2, f1, f2, null);
//				Comparison comp = ComparisonUtil.compareRatiosUsingSnpsFloat(ratios1, ratios2, f1, f2);
//				donorSB.append(comp.toString()).append("\n");
//				allComparisons.add(comp);
//			}
//			
//			// can now remove f1 from the cache as it is no longer required
////			Map<ChrPosition, float[]> m = cache.remove(f1);
//			TIntShortHashMap m = cache.remove(f1);
//			m.clear();
//			m = null;
//		}
		
		for (Comparison comp : allComparisons) {
			if (comp.getScore() > cutoff) {
				suspiciousResults.add(donor + "\t" + comp.toSummaryString());
			}
		}
		
		
		// flush out last donor details
		logger.info(donorSB.toString());
		
		logger.info("");
		if (suspiciousResults.isEmpty()) {
			logger.info("No suspicious results found");
		} else {
			logger.info("Suspicious results SUMMARY:");
			for (String s : suspiciousResults) logger.info(s);
		}
		
		if (outputXml != null)
			writeXmlOutput();
		
		return exitStatus;
	}
	
	private void performComparisons(List<File> files) {
		int size = files.size();
		AbstractQueue<Integer> queue =  new ConcurrentLinkedQueue<>();
		for (int i = 0 ; i < size -1 ; i++) {
			queue.add(i);
		}
		
		ExecutorService service = Executors.newFixedThreadPool(nThreads);		
		for (int i = 0 ; i < nThreads; i++) {
			service.execute(() -> {
					while (true) {
						Integer in = queue.poll();
						if (null == in) break;
				
						logger.info("performing comparison for : " + in.intValue());
						
						File f1 = files.get(in);
						TIntShortHashMap r1 = cache.get(f1);
						
						for (int j = in.intValue() + 1 ; j < size - 1 ; j ++ ) {
							File f2 = files.get(j);
							TIntShortHashMap r2 =  cache.get(f2);
							Comparison c = ComparisonUtil.compareRatiosUsingSnpsFloat(r1, r2, f1, f2);
							allComparisons.add(c);
						}
					}
				});
		}
		service.shutdown();
		try {
			if ( ! service.awaitTermination(100, TimeUnit.HOURS)) {
				logger.info("Timed out getting data from threads");
			}
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
	
	private void populateCache(List<File> files) {
		AbstractQueue<File> queue =  new ConcurrentLinkedQueue<>(files);
		ExecutorService service = Executors.newFixedThreadPool(nThreads);		
		for (int i = 0 ; i < nThreads; i++) {
			service.execute(() -> {
					while (true) {
						File f = queue.poll();
						if (null == f) break;
				
						logger.info("loading data from: " + f.getAbsolutePath());
						
						TIntShortHashMap genotypes = null;
						try {
							genotypes = SignatureUtil.loadSignatureRatiosFloatGenotype(f, minimumCoverage);
						} catch (Exception e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						TIntShortHashMap prevGenotypes = cache.putIfAbsent(f, genotypes);
						if (null != prevGenotypes) {
							logger.warn("already genotypes associated with file: " + f.getAbsolutePath());
						}
					}
				});
		}
		service.shutdown();
		try {
			if ( ! service.awaitTermination(100, TimeUnit.HOURS)) {
				logger.info("Timed out getting data from threads");
			}
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
//	TIntShortHashMap getSignatureData(File f) throws Exception {
//		// check map to see if this data has already been loaded
//		// if not - load
//		TIntShortHashMap result = cache.get(f);
//		if (result == null) {
//			result = SignatureUtil.loadSignatureRatiosFloatGenotype(f, minimumCoverage);
//			
//			if (result.size() < 1000) {
//				logger.warn("low coverage (" + result.size() + ") for file " + f.getAbsolutePath());
//			}
//			
////			if (cache.size() < cacheSize) {
//				cache.put(f, result);
////			}
//			fileIdsAndCounts.get(f)[1] = result.size();
//			/*
//			 * average coverage
//			 */
//			//TODO put this back in
////			IntSummaryStatistics iss = result.values().stream()
////				.mapToInt(array -> (int) array[4])
////				.summaryStatistics();
////			fileIdsAndCounts.get(f)[2] = (int) iss.getAverage();
//		}
//		return result;
//	}
//	Map<ChrPosition, float[]> getSignatureData(File f) throws Exception {
//		// check map to see if this data has already been loaded
//		// if not - load
//		Map<ChrPosition, float[]> result = cache.get(f);
//		if (result == null) {
//			result = SignatureUtil.loadSignatureRatiosFloat(f, minimumCoverage);
//			
//			if (result.size() < 1000) {
//				logger.warn("low coverage (" + result.size() + ") for file " + f.getAbsolutePath());
//			}
//			
//			if (cache.size() < cacheSize) {
//				cache.put(f, result);
//			}
//			fileIdsAndCounts.get(f)[1] = result.size();
//			/*
//			 * average coverage
//			 */
//			IntSummaryStatistics iss = result.values().stream()
//					.mapToInt(array -> (int) array[4])
//					.summaryStatistics();
//			fileIdsAndCounts.get(f)[2] = (int) iss.getAverage();
//		}
//		return result;
//	}
	
	private void writeXmlOutput() throws ParserConfigurationException, TransformerException {
		DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
		DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
		
		// root elements
		Document doc = docBuilder.newDocument();
		Element rootElement = doc.createElement("qsignature");
		doc.appendChild(rootElement);
		
		// list files
		Element filesE = doc.createElement("files");
		rootElement.appendChild(filesE);
		
		// write output xml file
		// do it to console first...
		List<File> keys = new ArrayList<>( fileIdsAndCounts.keySet());
		Collections.sort(keys, FileUtils.FILE_COMPARATOR);
		for (File f  : keys) {
			int[] value = fileIdsAndCounts.get(f);
			
			Element fileE = doc.createElement("file");
			fileE.setAttribute("id", value[0] + "");
			fileE.setAttribute("name", f.getAbsolutePath());
			fileE.setAttribute("coverage", value[1] + "");
			fileE.setAttribute("average_coverage_at_positions", value[2] + "");
			filesE.appendChild(fileE);
		}
		
		// list files
		Element compsE = doc.createElement("comparisons");
		rootElement.appendChild(compsE);
		for (Comparison comp : allComparisons) {
			int id1 = fileIdsAndCounts.get(comp.getMain())[0];
			int id2 = fileIdsAndCounts.get(comp.getTest())[0];
			
			Element compE = doc.createElement("comparison");
			compE.setAttribute("file1", id1 + "");
			compE.setAttribute("file2", id2 + "");
			compE.setAttribute("score", comp.getScore() + "");
			compE.setAttribute("overlap", comp.getOverlapCoverage() + "");
			compE.setAttribute("calcs", comp.getNumberOfCalculations() + "");
			compE.setAttribute("f1AveCovAtOverlaps", comp.getMainAveCovAtOverlaps() + "");
			compE.setAttribute("f2AveCovAtOverlaps", comp.getTestAveCovAtOverlaps() + "");
			compsE.appendChild(compE);
		}
		
		// write it out
		TransformerFactory transformerFactory = TransformerFactory.newInstance();
		Transformer transformer = transformerFactory.newTransformer();
		DOMSource source = new DOMSource(doc);
		StreamResult result = new StreamResult(new File(outputXml));
		transformer.setOutputProperty(OutputKeys.INDENT, "yes");
		transformer.transform(source, result);
	}
	
	private void addFilesToMap(List<File> orderedFiles) {
		int id = 1;
		for (File f : orderedFiles) {
			fileIdsAndCounts.put(f, new int[]{id++, -1, -1});
		}
	}

	
	public static void main(String[] args) throws Exception {
		LoadReferencedClasses.loadClasses(SignatureCompareRelatedSimpleGenotypeMT.class);
		
		SignatureCompareRelatedSimpleGenotypeMT sp = new SignatureCompareRelatedSimpleGenotypeMT();
		int exitStatus = 0;
		try {
			exitStatus = sp.setup(args);
		} catch (Exception e) {
			exitStatus = 2;
			if (null != logger)
				logger.error("Exception caught whilst running SignatureCompareRelatedSimple:", e);
			else {
				System.err.println("Exception caught whilst running SignatureCompareRelatedSimple: " + e.getMessage());
				System.err.println(Messages.USAGE);
			}
		}
		
		if (null != logger)
			logger.logFinalExecutionStats(exitStatus);
		
		System.exit(exitStatus);
	}
	
	protected int setup(String args[]) throws Exception{
		int returnStatus = 1;
		if (null == args || args.length == 0) {
			System.err.println(Messages.USAGE);
			System.exit(1);
		}
		Options options = new Options(args);

		if (options.hasHelpOption()) {
			System.err.println(Messages.USAGE);
			options.displayHelp();
			returnStatus = 0;
		} else if (options.hasVersionOption()) {
			System.err.println(Messages.getVersionMessage());
			returnStatus = 0;
		} else if ( ! options.hasLogOption()) {
			System.err.println(Messages.USAGE);
		} else {
			// configure logging
			logFile = options.getLog();
			logger = QLoggerFactory.getLogger(SignatureCompareRelatedSimpleGenotypeMT.class, logFile, options.getLogLevel());
			
			
			String [] cmdLineOutputFiles = options.getOutputFileNames();
			if (null != cmdLineOutputFiles && cmdLineOutputFiles.length > 0)
				outputXml = cmdLineOutputFiles[0];
			
			String[] paths = options.getDirNames(); 
			if (null != paths && paths.length > 0) {
				this.paths = paths;
			}
			if (null == paths) throw new QSignatureException("MISSING_DIRECTORY_OPTION");
			
			if (options.hasCutoff())
				cutoff = options.getCutoff();
			
			if (options.hasMinCoverage()) {
				minimumCoverage = options.getMinCoverage();
			}
			logger.tool("Setting minumim coverage to: " + minimumCoverage);
			
			additionalSearchStrings = options.getAdditionalSearchString();
			logger.tool("Setting additionalSearchStrings to: " + Arrays.deepToString(additionalSearchStrings));
			
			if (options.hasExcludeVcfsFileOption())
				excludeVcfsFile = options.getExcludeVcfsFile();
			
			logger.logInitialExecutionStats("SignatureCompareRelatedSimple", SignatureCompareRelatedSimpleGenotypeMT.class.getPackage().getImplementationVersion(), args);
			
			return engage();
		}
		return returnStatus;
	}

}