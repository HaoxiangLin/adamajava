/**
 * © Copyright The University of Queensland 2010-2014.  This code is released under the terms outlined in the included LICENSE file.
 */
package org.qcmg.qbamannotate;

import static java.util.Arrays.asList;

import java.util.List;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

public final class Options {
	private static final String HELP_DESCRIPTION = Messages
			.getMessage("HELP_OPTION_DESCRIPTION");
	private static final String VERSION_DESCRIPTION = Messages
			.getMessage("VERSION_OPTION_DESCRIPTION");
	private static final String SELF_DESCRIPTION = Messages
			.getMessage("SELF_OPTION_DESCRIPTION");
	private static final String TYPE_DESCRIPTION = Messages
			.getMessage("TYPE_OPTION_DESCRIPTION");
	private static final String LSIZE_DESCRIPTION = Messages
			.getMessage("LSIZE_OPTION_DESCRIPTION");
	private static final String USIZE_DESCRIPTION = Messages
			.getMessage("USIZE_OPTION_DESCRIPTION");
	private static final String STACKTRACE_OPTION_DESCRIPTION = Messages
			.getMessage("STACKTRACE_OPTION_DESCRIPTION");
	private static final String LOG_OPTION_DESCRIPTION = Messages
			.getMessage("LOG_OPTION_DESCRIPTION");
	private static final String LOG_LEVEL_OPTION_DESCRIPTION = Messages
			.getMessage("LOG_LEVEL_OPTION_DESCRIPTION");

	private final OptionParser parser = new OptionParser();
	private final OptionSet options;
	private final boolean hasInsufficientArguments;
	private final boolean hasTooManyArguments;
	private final String outputFileName;
	private final String inputBAMFileName;
	private final String inputMAFileNameA;
	private final String inputMAFileNameB;
	private final Integer[] lowerISizeLimits;
	private final Integer[] upperISizeLimits;
	private final String[] types;
	private final String log;
	private final String logLevel;

	@SuppressWarnings("unchecked")
	public Options(final String[] args) throws Exception {
		parser.acceptsAll(asList("h", "help"), HELP_DESCRIPTION);
		parser.acceptsAll(asList("v", "V", "version"), VERSION_DESCRIPTION);
		parser.acceptsAll(asList("r", "reportbams"), SELF_DESCRIPTION);
		parser.acceptsAll(asList("s", "self"), SELF_DESCRIPTION);
		parser.acceptsAll(asList("t", "type"), TYPE_DESCRIPTION)
				.withRequiredArg().ofType(String.class).describedAs("run type");
		parser.acceptsAll(asList("l", "lsize"), LSIZE_DESCRIPTION)
				.withRequiredArg().ofType(Integer.class).describedAs(
						"lower isize");
		parser.acceptsAll(asList("u", "usize"), USIZE_DESCRIPTION)
				.withRequiredArg().ofType(Integer.class).describedAs(
						"upper isize");
		parser.acceptsAll(asList("T", "stack-trace"),
				STACKTRACE_OPTION_DESCRIPTION);
		parser.accepts("log", LOG_OPTION_DESCRIPTION).withRequiredArg().ofType(
				String.class);
		parser.accepts("loglevel", LOG_LEVEL_OPTION_DESCRIPTION)
				.withRequiredArg().ofType(String.class);
		parser.posixlyCorrect(true);
		options = parser.parse(args);
		List<String> nonoptions = options.nonOptionArguments();
		if (!hasSelfOption()) {
			if (3 > nonoptions.size()) {
				hasInsufficientArguments = true;
			} else {
				hasInsufficientArguments = false;
			}
			if (4 < nonoptions.size()) {
				hasTooManyArguments = true;
			} else {
				hasTooManyArguments = false;
			}
			if (0 < nonoptions.size()) {
				outputFileName = nonoptions.get(0).trim();
			} else {
				outputFileName = null;
			}
			if (1 < nonoptions.size()) {
				inputBAMFileName = nonoptions.get(1).trim();
			} else {
				inputBAMFileName = null;
			}
			if (2 < nonoptions.size()) {
				inputMAFileNameA = nonoptions.get(2).trim();
			} else {
				inputMAFileNameA = null;
			}
			if (3 < nonoptions.size()) {
				inputMAFileNameB = nonoptions.get(3).trim();
			} else {
				inputMAFileNameB = null;
			}
		} else {
			if (2 > nonoptions.size()) {
				hasInsufficientArguments = true;
			} else {
				hasInsufficientArguments = false;
			}
			if (3 < nonoptions.size()) {
				hasTooManyArguments = true;
			} else {
				hasTooManyArguments = false;
			}
			if (0 < nonoptions.size()) {
				outputFileName = nonoptions.get(0).trim();
			} else {
				outputFileName = null;
			}
			if (1 < nonoptions.size()) {
				inputBAMFileName = nonoptions.get(1).trim();
			} else {
				inputBAMFileName = null;
			}
			inputMAFileNameA = null;
			inputMAFileNameB = null;
		}

		List typesList = options.valuesOf("t");
		types = new String[typesList.size()];
		typesList.toArray(types);

		List lSizeList = options.valuesOf("l");
		lowerISizeLimits = new Integer[lSizeList.size()];
		lSizeList.toArray(lowerISizeLimits);

		List uSizeList = options.valuesOf("u");
		upperISizeLimits = new Integer[uSizeList.size()];
		uSizeList.toArray(upperISizeLimits);

		log = (String) options.valueOf("log");
		logLevel = (String) options.valueOf("loglevel");
	}

	public boolean hasStatsDirOption() {
		return options.has("D") || options.has("stats-dir");
	}

	public boolean hasStatsBasenameOption() {
		return options.has("B") || options.has("stats-basename");
	}

	public boolean hasStatsMatchBasename() {
		return options.has("M") || options.has("stats-match-basename");
	}

	String getLog() {
		return log;
	}

	String getLogLevel() {
		return logLevel;
	}

	public boolean hasMalformedOption() {
		return options.has("m") || options.has("malformed");
	}

	public boolean hasPerFeatureOption() {
		return options.has("per-feature");
	}

	public boolean hasStackTraceOption() {
		return options.has("T") || options.has("stack-trace");
	}

	public boolean hasVersionOption() {
		return options.has("v") || options.has("V") || options.has("version");
	}

	public boolean hasHelpOption() {
		return options.has("h") || options.has("help");
	}

	public boolean hasSelfOption() {
		return options.has("s") || options.has("self");
	}

	public boolean hasReportBamsOption() {
		return options.has("r") || options.has("reportbams");
	}

	public boolean hasTypeOption() {
		return options.has("t") || options.has("type");
	}

	public boolean hasLSizeOption() {
		return options.has("l") || options.has("lsize");
	}

	public boolean hasUSizeOption() {
		return options.has("u") || options.has("usize");
	}

	public boolean hasInsufficientArguments() {
		return hasInsufficientArguments;
	}

	public boolean hasTooManyArguments() {
		return hasTooManyArguments;
	}

	public String getOutputFileName() {
		return outputFileName;
	}

	public String getInputBAMFileName() {
		return inputBAMFileName;
	}

	public boolean hasInputMAFileNameA() {
		return null != inputMAFileNameA;
	}

	public boolean hasInputMAFileNameB() {
		return null != inputMAFileNameB;
	}

	boolean hasNumberRecordsOption() {
		return options.has("n") || options.has("number");
	}

	public String getInputMAFileNameA() {
		return inputMAFileNameA;
	}

	public String getInputMAFileNameB() {
		return inputMAFileNameB;
	}

	public void displayHelp() throws Exception {
		parser.printHelpOn(System.out);
	}

	public Integer[] getLowerISizeLimit() {
		return lowerISizeLimits;
	}

	public Integer[] getUpperISizeLimit() {
		return upperISizeLimits;
	}

	public boolean hasTooManyTypeArguments() {
		return types.length > 1;
	}

	public boolean hasTypeArgument() {
		return types.length > 0;
	}

	public Object[] getTypes() {
		return types;
	}

	public boolean hasLSizeLimit() {
		return lowerISizeLimits.length > 0;
	}

	public boolean hasTooManyLSizeLimits() {
		return lowerISizeLimits.length > 1;
	}

	public Integer[] getLowerISizeLimits() {
		return lowerISizeLimits;
	}

	public boolean hasUSizeLimit() {
		return upperISizeLimits.length > 0;
	}

	public boolean hasTooManyUSizeLimits() {
		return upperISizeLimits.length > 1;
	}

	public Integer[] getUpperISizeLimits() {
		return upperISizeLimits;
	}

	public boolean hasInputMAFileNames() {
		return null != inputMAFileNameA || null != inputMAFileNameB;
	}

	public void detectBadOptions() throws Exception {
		if (hasInsufficientArguments()) {
			throw new BamAnnotateException("INSUFFICIENT_ARGUMENTS");
		}
		if (hasTooManyArguments()) {
			throw new BamAnnotateException("TOO_MANY_ARGUMENTS");
		}
		if (!hasLogOption()) {
			throw new Exception(
					"A log filename must be specified (using the --log option)");
		}
	}

	public boolean hasLogOption() {
		return options.has("log");
	}

	public boolean hasLogLevelOption() {
		return options.has("loglevel");
	}
}