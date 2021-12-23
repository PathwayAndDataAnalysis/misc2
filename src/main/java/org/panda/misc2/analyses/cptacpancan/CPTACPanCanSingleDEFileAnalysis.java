package org.panda.misc2.analyses.cptacpancan;

import org.panda.utility.*;
import org.panda.utility.statistics.UniquePrinter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.*;
import java.util.stream.Collectors;

public class CPTACPanCanSingleDEFileAnalysis
{
	public static final String DATA_BASE = "/home/ozgunbabur/Data/CPTAC-PanCan/";
	public static final String OUT_BASE = "/home/ozgunbabur/Analyses/CPTAC-PanCan/";
	public static final String MAP_FILE = DATA_BASE + "var_map_full.tsv";

	public static final String DELIM = "\t";

	static List<String> KEGG_NAMES = Arrays.asList("KEGG_FATTY_ACID_METABOLISM", "KEGG_GLYCOLYSIS_GLUCONEOGENESIS", "KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_FOCAL_ADHESION");

	public static void main(String[] args) throws IOException
	{
		convertDiffExpData(DATA_BASE + "DDR_res.tsv", OUT_BASE + "DDR_res.tsv");
		prepareAnalysisFoldersDiffExp(OUT_BASE + "DDR_res.tsv", OUT_BASE + "DDR_res/");

		convertDiffExpData(DATA_BASE + "DDRSub_res.tsv", OUT_BASE + "DDRSub_res.tsv");
		prepareAnalysisFoldersDiffExp(OUT_BASE + "DDRSub_res.tsv", OUT_BASE + "DDRSub_res/");
	}

	public static void process(String inFile, String inBase, String outBase) throws IOException
	{
		String inDir = inFile.substring(0, inFile.lastIndexOf(File.separator) + 1);
		String name = inFile.substring(inFile.lastIndexOf(File.separator) + 1, inFile.lastIndexOf("."));

		if (!inDir.startsWith(inBase)) throw new RuntimeException("The file has to be nested under the given inBase");

		String outDir = inDir.replaceFirst(inBase, outBase) + File.separator + name;
		String outDataFile = outDir + "/data.txt";

		convertDiffExpData(inFile, outDataFile);
		prepareAnalysisFoldersDiffExp(outDataFile, outDir);
	}

	private static void convertDiffExpData(String inFile, String outFile) throws IOException
	{
		FileUtil.mkdirsOfFilePath(outFile);

//		Map<String, String> idToGene = readIDToGeneMapping();

		String[] header = FileUtil.lines(inFile).findFirst().get().split("\t");
		int idInd = ArrayUtil.indexOf(header, "index");
		int geneInd = ArrayUtil.indexOf(header, "gene_name");
		int featInd = ArrayUtil.indexOf(header, "feature");
		int pInd = ArrayUtil.indexOf(header, "adj.P.Val");
		int fcInd = ArrayUtil.indexOf(header, "logFC");
		int clusterInd = ArrayUtil.indexOf(header, "id");

		Map<String, Map<String, Double>> idToVals = new HashMap<>();
		Map<String, String[]> idToProps = new HashMap<>();

		Set<String> mem = new HashSet<>();

		UniquePrinter up = new UniquePrinter();
		FileUtil.lines(inFile).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			if (t[pInd].isEmpty()) return;

			String id = t[idInd];

			String[] props = idToProps.get(id);

			if (props == null)
			{
//				String sym = idToGene.get(id);
				String sym = t[geneInd];
//				if (HGNC.get().getSymbol(sym) == null) return;
				String mod;
				switch (t[featInd])
				{
					case "phosphoproteome":
					case "phosphoproteome_res":
						mod = "P";
						break;
					case "acetylome":
					case "acetylome_res":
						mod = "A";
						break;
					case "transcriptome":
						mod = "R";
						break;
					case "proteome":
						mod = "G";
						break;
					default:
						up.print("Unknown feature: ", t[featInd]);
						mod = "";
				}

				if (mod.equals("P") || mod.equals("A"))
				{
					List<String> siteList = getSitesFromID(id);

					String site = CollectionUtil.merge(siteList, "|");

					String cpID = sym + "-" + mod + "-" + site.replaceAll("\\|", "-");

					while (mem.contains(cpID))
					{
						String ss = cpID.substring(cpID.lastIndexOf("-") + 1);
						cpID = Character.isDigit(ss.charAt(0)) ? cpID.substring(0, cpID.lastIndexOf("-") + 1) + (Integer.valueOf(ss) + 1) : cpID + "-2";
					}
					mem.add(cpID);

					props = new String[]{cpID, sym, site, mod, ""};
				}
				else
				{
					props = new String[]{sym + "-" + mod, sym, "", mod, ""};
				}
				idToProps.put(id, props);
			}

			double p = Double.valueOf(t[pInd]);
			if (t[fcInd].startsWith("-")) p = -p;

			if (!idToVals.containsKey(id)) idToVals.put(id, new HashMap<>());
			String cluster = clusterInd < 0 ? "single-case" : t[clusterInd];
			idToVals.get(id).put(cluster, p);
		});

		List<String> clusters = idToVals.keySet().stream().map(idToVals::get).map(Map::keySet).flatMap(Collection::stream).distinct().collect(Collectors.toList());
		clusters.sort(String::compareTo);
//		clusters.sort(Comparator.comparing(Double::valueOf));

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("ID\tSymbols\tSites\tFeature\tEffect");
		clusters.forEach(c -> FileUtil.tab_write(c, writer));
		idToProps.forEach((id, props) ->
		{
			FileUtil.lnwrite(ArrayUtil.merge("\t", props), writer);
			Map<String, Double> valMap = idToVals.get(id);
			clusters.forEach(c -> FileUtil.tab_write(valMap.getOrDefault(c, Double.NaN), writer));
		});
		writer.close();
	}

	private static void prepareAnalysisFoldersDiffExp(String dataFile, String rootDir) throws IOException
	{
		FileUtil.mkdirs(rootDir);
		String dataName = dataFile.substring(dataFile.lastIndexOf("/") + 1);

		List<String> clusters = Arrays.asList(FileUtil.lines(dataFile).findFirst().get().split("\t"));
		clusters = clusters.subList(5, clusters.size());

		for (int i = 0; i < clusters.size(); i++)
		{
			String c1 = clusters.get(i);
			String dir = rootDir + File.separator + c1 + File.separator;
			FileUtil.mkdirs(dir);

			String valPointers = "value-column = " + c1;

			String innerDir = dir + "sitespec/";
			new File(innerDir).mkdirs();
			BufferedWriter writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
			writer.write(getSiteSpecParams(dataName) + valPointers);
			writer.close();
			innerDir = dir + "globprot/";
			new File(innerDir).mkdirs();
			writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
			writer.write(getGlobProtParams(dataName) + valPointers);
			writer.close();
			innerDir = dir + "rnaseq/";
			new File(innerDir).mkdirs();
			writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
			writer.write(getRNAParams(dataName) + valPointers);
			writer.close();
			innerDir = dir + "all/";
			new File(innerDir).mkdirs();
			writer = new BufferedWriter(new FileWriter(innerDir + "parameters.txt"));
			writer.write(getAllParams(dataName) + valPointers);
			writer.close();
		}
	}


	private static Map<String, String> readIDToGeneMapping() throws IOException
	{
		return Files.lines(Paths.get(MAP_FILE)).skip(1)
			.map(l -> l.split("\t"))
//			.peek(t -> t[0] = t[0].startsWith("ENSG") ? t[0].substring(0, t[0].indexOf(".")) : t[0])
			.collect(Collectors.toMap(t -> t[0], t -> t[2]));//, (o, o2) -> o));
	}

	private static List<String> getSitesFromID(String id)
	{
		id = id.split("_")[2];
		id = id.substring(0, id.length() - 1);
		return Arrays.asList(id.split("s|t|y|k"));
	}

	public static String getParamStart(String dataName)
	{
		return "proteomics-values-file = ../../" + dataName + "\n" +
			"id-column = ID\n" +
			"symbols-column = Symbols\n" +
			"sites-column = Sites\n" +
			"feature-column = Feature\n" +
			"effect-column = Effect\n" +
			"\n" +
			"value-transformation = signed-p-values\n" +
			"threshold-for-data-significance = 0.1 protein\n" +
			"threshold-for-data-significance = 0.1 phosphoprotein\n" +
			"threshold-for-data-significance = 0.1 acetylprotein\n" +
			"threshold-for-data-significance = 0.1 rna\n" +
			"\n" +
			"color-saturation-value = 15\n" +
			"\n" +
			"calculate-network-significance = true\n" +
			"permutations-for-significance = 10000\n" +
			"fdr-threshold-for-network-significance = 0.1\n" +
			"use-network-significance-for-causal-reasoning = true\n" +
			"\n";
	}

	public static String getSiteSpecParams(String dataName)
	{
		return getParamStart(dataName) +
			"relation-filter-type = site-specific-only\n" +
			"\n";
	}

	public static String getGlobProtParams(String dataName)
	{
		return getParamStart(dataName) +
			"relation-filter-type = expression-only\n" +
			"data-type-for-expressional-targets = protein\n" +
			"\n";
	}

	public static String getAllParams(String dataName)
	{
		return getParamStart(dataName) +
			"data-type-for-expressional-targets = rna\n" +
			"\n";
	}

	public static String getRNAParams(String dataName)
	{
		return getAllParams(dataName) +
			"relation-filter-type = expression-only\n" +
			"\n";
	}
}
