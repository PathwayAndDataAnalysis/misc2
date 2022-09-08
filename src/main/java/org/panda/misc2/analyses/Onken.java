package org.panda.misc2.analyses;

import org.panda.misc2.KinaseEnrichment;
import org.panda.misc2.causalpath.CausalPathSubnetwork;
import org.panda.resource.GO;
import org.panda.resource.HGNC;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.StreamDirection;
import org.panda.utility.UniqueMaker;
import org.panda.utility.statistics.RankedListSignedGroupedDifferentialEnrichment;
import org.panda.utility.statistics.RankedListSignedGroupedEnrichment;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

public class Onken {
//		public static final String DATA_DIR = "/home/ozgunbabur/Data/Onken/";
	public static final String DATA_DIR = "/Users/ozgun/Data/Onken/";
//		public static final String OUT_DIR = "/home/ozgunbabur/Analyses/Onken/";
	public static final String OUT_DIR = "/Users/ozgun/Analyses/Onken/";

	public static final String[] CELL_LINES = new String[]{"MP41", "MP46", "OCM"};

	public static void main(String[] args) throws IOException {
		convertData();
//		predictKinaseActivity();
//		generateSubgraphs();
	}

	private static void generateSubgraphs() throws IOException
	{
		Map<String, Set<String>> goiMap = getGOIMap();

		CausalPathSubnetwork.writeGOINeighForCompBasedRecursively(OUT_DIR, goiMap, StreamDirection.BOTHSTREAM);
	}

	public static void predictKinaseActivity() throws IOException {
		Map<String, List<String>> rankedLists = new HashMap<>();

		for (String cellLine : CELL_LINES) {
			List<String> list = KinaseEnrichment.readRankedIDsFromCPFile(OUT_DIR + cellLine + "/data.tsv");
			rankedLists.put(cellLine, list);
		}

		Map<String, Map<String, Map<String, Boolean>>> rawPriors = KinaseEnrichment.readPriors();

		for (String cl1 : rankedLists.keySet()) {

			String dataFile1 = OUT_DIR + cl1 + "/data.tsv";
			Map<String, Map<String, Set<String>>> geneToSiteToID = KinaseEnrichment.readDataMappingFromCPFile(dataFile1);
			Map<String, Set<Map<String, Boolean>>> priors = KinaseEnrichment.convertPriors(geneToSiteToID, rawPriors, rankedLists.get(cl1));

			RankedListSignedGroupedEnrichment.reportEnrichment(rankedLists.get(cl1), priors, 1000000, dataFile1.substring(0, dataFile1.lastIndexOf("/")) + "/kinase-enrichment.tsv");

			for (String cl2 : rankedLists.keySet()) {
				if (!cl2.equals(cl1)) {
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
			convertOneCellLine(cellLine);
		}
	}

	public static void convertOneCellLine(String cellLine) throws IOException
	{
		String dir = OUT_DIR + cellLine + "/";
		FileUtil.mkdirs(dir);
		BufferedWriter writer = FileUtil.newBufferedWriter(dir + "data.tsv");
		writer.write("ID\tSymbols\tSites\tFeature\tEffect\tSignedP");
		writer.write(convertPhospho(DATA_DIR + cellLine + ".csv"));
		writer.write(convertTotProt(DATA_DIR + cellLine + "-total-protein.csv"));
		if (!cellLine.equals("OCM"))
			writer.write(convertRNA(DATA_DIR + cellLine + "-rna-seq.csv"));
		writer.close();
	}

	public static String convertPhospho(String inFile) {
		StringBuilder sb = new StringBuilder();

		String[] header = FileUtil.readHeader(inFile);
		int symInd = ArrayUtil.indexOf(header, "Gene.names");
		int altSymInd = ArrayUtil.indexOf(header, "GENE");
		int siteInd = ArrayUtil.indexOf(header, "MOD_RSD");
		int fcInd = ArrayUtil.indexOf(header, "log2FC");
		int pInd = ArrayUtil.indexOf(header, "pvalue");

		UniqueMaker um = new UniqueMaker();

		FileUtil.linesTabbedSkip1(inFile).forEach(t ->
		{
			String sym = t[symInd];
			String site = t[siteInd];

			String id = um.get(sym + "-" + site);

			double signedP = Double.parseDouble(t[pInd]);
			if (t[fcInd].startsWith("-")) signedP = -signedP;

			sb.append("\n").append(id).append("\t").append(sym).append("\t").append(site).append("\t").append("P").append("\t\t").append(signedP);

			if (t.length > altSymInd) {
				String altSym = t[altSymInd];
				if (!altSym.isEmpty() && !altSym.equals(sym)) {
					System.out.println(sym + ":" + HGNC.get().getSymbol(sym) + "\t" + altSym + ":" + HGNC.get().getSymbol(altSym));
				}
			}
		});
		return sb.toString();
	}

	public static String convertTotProt(String inFile) {
		StringBuilder sb = new StringBuilder();

		String[] header = FileUtil.readHeader(inFile);
		int symInd = ArrayUtil.indexOf(header, "Gene.names");
		int fcInd = ArrayUtil.indexOf(header, "log2FC");
		int pInd = ArrayUtil.indexOf(header, "pvalue");

		UniqueMaker um = new UniqueMaker();

		FileUtil.linesTabbedSkip1(inFile).forEach(t ->
		{
			String sym = t[symInd];

			String id = um.get(sym);

			double signedP = Double.parseDouble(t[pInd]);
			if (t[fcInd].startsWith("-")) signedP = -signedP;

			sb.append("\n").append(id).append("\t").append(sym).append("\t\tG\t\t").append(signedP);
		});
		return sb.toString();
	}

	public static String convertRNA(String inFile) {
		StringBuilder sb = new StringBuilder();

		String[] header = FileUtil.readHeader(inFile);
		int symInd = ArrayUtil.indexOf(header, "external_gene_name");
		int typeInd = ArrayUtil.indexOf(header, "external_gene_source");
		int fcInd = ArrayUtil.indexOf(header, "logFC");
		int pInd = ArrayUtil.indexOf(header, "pval");

		UniqueMaker um = new UniqueMaker();

		FileUtil.linesTabbedSkip1(inFile).forEach(t ->
		{
			String sym = t[symInd];

			if (!t[typeInd].equals("HGNC Symbol")) return;

			String id = um.get(sym + "-rna");

			double signedP = Double.parseDouble(t[pInd]);
			if (t[fcInd].startsWith("-")) signedP = -signedP;

			sb.append("\n").append(id).append("\t").append(sym).append("\t\tR\t\t").append(signedP);
		});
		return sb.toString();
	}

	private static Map<String, Set<String>> getGOIMap()
	{
		Map<String, Set<String>> goiMap = new HashMap<>();
		goiMap.put("CDK1-2", new HashSet<>(Arrays.asList("CDK1", "CDK2")));
		goiMap.put("RPS6KA1", new HashSet<>(List.of("RPS6KA1")));
		goiMap.put("PRKC", new HashSet<>(Arrays.asList("PRKCA", "PRKCD")));
		goiMap.put("AURKA", new HashSet<>(List.of("AURKA")));
		goiMap.put("MAPK1-3", new HashSet<>(Arrays.asList("MAPK1", "MAPK3")));
		goiMap.put("CSNK2A1", new HashSet<>(List.of("CSNK2A1")));
		goiMap.put("PRKACA", new HashSet<>(List.of("PRKACA")));
		goiMap.put("PRKA", new HashSet<>(Arrays.asList("PRKAB1", "PRKAG2")));
		goiMap.put("RPS6K", new HashSet<>(Arrays.asList("RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KA6", "RPS6KB1", "RPS6KB2")));
		goiMap.put("CDK7-12-13", new HashSet<>(Arrays.asList("CDK7", "CDK12", "CDK13")));

		goiMap.put("carbohydrate-metabolic-process", GO.get().getGenesOfTerm("GO:0005975"));
		goiMap.put("regulation-of-carbohydrate-metabolic-process", GO.get().getGenesOfTerm("GO:0006109"));
		goiMap.put("lipid-metabolic-process", GO.get().getGenesOfTerm("GO:0006629"));
		goiMap.put("regulation-of-lipid-metabolic-process", GO.get().getGenesOfTerm("GO:0019216"));
		goiMap.put("tricarboxylic-acid-cycle", GO.get().getGenesOfTerm("GO:0006099"));

		return goiMap;
	}
}