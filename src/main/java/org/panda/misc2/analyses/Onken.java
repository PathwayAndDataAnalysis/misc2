package org.panda.misc2.analyses;

import org.panda.misc2.KinaseEnrichment;
import org.panda.resource.HGNC;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.UniqueMaker;
import org.panda.utility.statistics.RankedListSignedGroupedDifferentialEnrichment;
import org.panda.utility.statistics.RankedListSignedGroupedEnrichment;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

public class Onken
{
	public static final String DATA_DIR = "/home/ozgunbabur/Data/Onken/";
	public static final String OUT_DIR = "/home/ozgunbabur/Analyses/Onken/";

	public static final String[] CELL_LINES = new String[]{"MP41", "MP46", "OCM"};

	public static void main(String[] args) throws IOException
	{
		convertData();
//		predictKinaseActivity();
	}

	public static void predictKinaseActivity() throws IOException
	{
		Map<String, List<String>> rankedLists = new HashMap<>();

		for (String cellLine : CELL_LINES)
		{
			List<String> list = KinaseEnrichment.readRankedIDsFromCPFile(OUT_DIR + cellLine + "/data.tsv");
			rankedLists.put(cellLine, list);
		}

		Map<String, Map<String, Map<String, Boolean>>> rawPriors = KinaseEnrichment.readPriors();

		for (String cl1 : rankedLists.keySet())
		{

			String dataFile1 = OUT_DIR + cl1 + "/data.tsv";
			Map<String, Map<String, Set<String>>> geneToSiteToID = KinaseEnrichment.readDataMappingFromCPFile(dataFile1);
			Map<String, Set<Map<String, Boolean>>> priors = KinaseEnrichment.convertPriors(geneToSiteToID, rawPriors, rankedLists.get(cl1));

			RankedListSignedGroupedEnrichment.reportEnrichment(rankedLists.get(cl1), priors, 1000000, dataFile1.substring(0, dataFile1.lastIndexOf("/")) + "/kinase-enrichment.tsv");

			for (String cl2 : rankedLists.keySet())
			{
				if (!cl2.equals(cl1))
				{
					System.out.println(cl1 + " -vs- " + cl2);
					List<String> rankedList1 = new ArrayList<>(rankedLists.get(cl1));
					List<String> rankedList2 = new ArrayList<>(rankedLists.get(cl2));
					System.out.println("rankedList1.size() = " + rankedList1.size());
					System.out.println("rankedList2.size() = " + rankedList2.size());

					rankedList1.retainAll(rankedList2);
					rankedList2.retainAll(rankedList1);
					System.out.println("common size        = " + rankedList2.size());

					geneToSiteToID = KinaseEnrichment.readDataMappingFromCPFile(dataFile1);
					priors = KinaseEnrichment.convertPriors(geneToSiteToID, rawPriors, rankedList1);

					String outFile = OUT_DIR + cl1 + "-vs-" + cl2 + "-differential-kinase-enrichment.tsv";

					RankedListSignedGroupedDifferentialEnrichment.reportEnrichment(rankedList1, rankedList2, priors, 1000000, outFile);
				}
			}
		}
	}


	public static void convertData() throws IOException
	{
		for (String cellLine : CELL_LINES)
		{
			String dir = OUT_DIR + cellLine + "/";
			FileUtil.mkdirs(dir);
			convertSinglePhospho(DATA_DIR + cellLine + ".csv", dir + "data.tsv");
		}
	}

	public static void convertSinglePhospho(String inFile, String outFile) throws IOException
	{
		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("ID\tSymbols\tSites\tFeature\tEffect\tSignedP");
		writer.write(convertPhospho(inFile));
		writer.close();
	}

	public static String convertPhospho(String inFile)
	{
		StringBuilder sb = new StringBuilder();

		String[] header = FileUtil.readHeader(inFile);
		int symInd = ArrayUtil.indexOf(header, "Gene.names");
		int altSymInd = ArrayUtil.indexOf(header, "GENE");
		int siteInd = ArrayUtil.indexOf(header, "MOD_RSD");
		int fcInd = ArrayUtil.indexOf(header, "log2FC");
		int pInd = ArrayUtil.indexOf(header, "adj.pvalue");

		UniqueMaker um = new UniqueMaker();

		FileUtil.linesTabbedSkip1(inFile).forEach(t ->
		{
			String sym = t[symInd];
			String site = t[siteInd];

			String id = um.get(sym + "-" + site);

			double signedP = Double.parseDouble(t[pInd]);
			if (t[fcInd].startsWith("-")) signedP = -signedP;

			sb.append("\n").append(id).append("\t").append(sym).append("\t").append(site).append("\t").append("P").append("\t\t").append(signedP);

			if (t.length > altSymInd)
			{
				String altSym = t[altSymInd];
				if (!altSym.isEmpty() && !altSym.equals(sym))
				{
					System.out.println(sym + ":" + HGNC.get().getSymbol(sym) + "\t" + altSym + ":" + HGNC.get().getSymbol(altSym));
				}
			}
		});
		return sb.toString();
	}
}
