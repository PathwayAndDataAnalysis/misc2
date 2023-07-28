package org.panda.misc2;

import org.panda.utility.FileUtil;
import org.panda.utility.statistics.RankedListSignedGroupedDifferentialEnrichment;
import org.panda.utility.statistics.RankedListSignedGroupedEnrichment;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class DiferentialKinaseEnrichment extends KinaseEnrichment
{
	public static void main(String[] args) throws IOException
	{
		analyzePARThrombin();
//		analyzeCond1Cond2();
	}

	public static void analyzePARThrombin() throws IOException
	{
		Map<String, Map<String, Map<String, Boolean>>> rawPriors = readPriors();

		String dir = "/home/ozgunbabur/Analyses/Aslan-Thrombin-PAR/";
		String dataFile = dir + "data.csv";
		List<String> rankedListControl = readRankedIDsFromCPFile(dataFile, "Resting-vs-PAR1");
		List<String> rankedListTest = readRankedIDsFromCPFile(dataFile, "Resting-vs-Thrombin");


		Map<String, Map<String, Set<String>>> geneToSiteToID = readDataMappingFromCPFile(dataFile);
		filterPriorsToExisting(rawPriors, dataFile);
		Map<String, Set<Map<String, Boolean>>> priors = convertPriors(geneToSiteToID, rawPriors, rankedListControl);

		String outFile = dir + "WithKinaseLibrary/diff-kinase-enrichment-Thrombin-vs-PAR1.tsv";
		RankedListSignedGroupedDifferentialEnrichment.reportEnrichment(rankedListTest, rankedListControl, priors, 1000000, outFile);
	}

	public static void analyzeCond1Cond2() throws IOException
	{
		Map<String, Map<String, Map<String, Boolean>>> rawPriors = readPriors();

		String dir = "/home/ozgunbabur/Data/Aslan/gpvi_mass_spec/";
		String dataFile1 = dir + "C1S.tsv";
		String dataFile2 = dir + "C2S.tsv";
		List<String> rankedListControl = readRankedIDsFromGP6Sheet(dataFile1);
		List<String> rankedListTest = readRankedIDsFromGP6Sheet(dataFile2);

		rankedListControl.retainAll(rankedListTest);
		rankedListTest.retainAll(rankedListControl);

		Map<String, Map<String, Set<String>>> geneToSiteToID = readDataMappingFromGP6Sheet(dataFile1);
//		filterPriorsToExisting(rawPriors, dataFile1, dataFile2);
		Map<String, Set<Map<String, Boolean>>> priors = convertPriors(geneToSiteToID, rawPriors, rankedListControl);

		String outFile = dir + "diff-kinase-enrichment-with-KL-cond2-vs-cond1.tsv";
		RankedListSignedGroupedDifferentialEnrichment.reportEnrichment(rankedListTest, rankedListControl, priors, 1000000, outFile);
	}

	public static Map<String, Map<String, Map<String, Boolean>>> filterPriorsToExisting(Map<String,
		Map<String, Map<String, Boolean>>> rawPriors, String cpDataFile1, String cpDataFile2)
	{
		Set<String> existing = FileUtil.getTermsInTabDelimitedColumn(cpDataFile1, 1, 1);
		existing.addAll(FileUtil.getTermsInTabDelimitedColumn(cpDataFile2, 1, 1));
		Set<String> remove = rawPriors.keySet().stream().filter(k -> !existing.contains(k)).collect(Collectors.toSet());
		remove.forEach(rawPriors::remove);
		return rawPriors;
	}

}
