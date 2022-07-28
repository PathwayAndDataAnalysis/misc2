package org.panda.misc2.analyses;

import org.panda.causalpath.network.GraphWriter;
import org.panda.causalpath.run.JasonizeResultGraphsRecursively;
import org.panda.resource.MGI;
import org.panda.resource.MSigDB;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;
import org.panda.utility.statistics.RankedListDifferentialEnrichment;
import org.panda.utility.statistics.RankedListEnrichment;
import org.panda.utility.statistics.RankedListSignedDifferentialEnrichment;
import org.panda.utility.statistics.RankedListSignedEnrichment;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

public class Josh
{
	public static final String DATA_DIR = "/home/ozgunbabur/Data/Josh/cluster_list/";
	public static final String OUT_DIR = "/home/ozgunbabur/Analyses/Josh/";
	static int RAND_ITER = 1000000;

	public static void main(String[] args) throws IOException
	{
		String[] filename = new String[]{"4nqo_epithelial.xls", "M-TPC.xls", "T-TPC.xls", "tumor.xls", "vehicle_epithelial.xls"};
		Map<String, Map<String, Boolean>> priors = readTFTargetPriors();

		Map<String, List<String>> rankedLists = new HashMap<>();
		Map<String, List<String>> signedRankedLists = new HashMap<>();

		for (String file : filename)
		{
			System.out.println("file = " + file);
			String fileNoExt = file.substring(0, file.lastIndexOf("."));
			FileUtil.mkdirs(OUT_DIR + fileNoExt);

			List<String> mouseList = readRankedList(DATA_DIR + file);
			List<String> humanList = convertToHuman(mouseList);
			rankedLists.put(fileNoExt, humanList);

			System.out.println("rankedList = " + mouseList.size());
			System.out.println("humanList  = " + humanList.size());

			// for sanity check
//			Collections.shuffle(humanList);

//			Map<String, Set<String>> geneSets = loadReactomeGeneSets();
//			reportEnrichments(humanList, cropGeneSetsToSpecificList(geneSets, humanList), OUT_DIR + fileNoExt + "/" + fileNoExt + ".reactome.xls");
//			geneSets = loadHallmarkGeneSets();
//			reportEnrichments(humanList, cropGeneSetsToSpecificList(geneSets, humanList), OUT_DIR + fileNoExt + "/" +fileNoExt + ".hallmark.xls");

			mouseList = readSignedRankedList(DATA_DIR + file);
			humanList = convertToHuman(mouseList);
			signedRankedLists.put(fileNoExt, humanList);
			System.out.println("rankedList = " + mouseList.size());
			System.out.println("humanList  = " + humanList.size());

			Map<String, Map<String, Boolean>> croppedPriors = cropPriorsToSpecificList(priors, humanList);
			String tfResultFile = OUT_DIR + fileNoExt + "/" + fileNoExt + ".TF-activity.xls";
//			reportTFTargetSignedEnricments(humanList, croppedPriors, tfResultFile);
//			generateNetworkResults(DATA_DIR + file, tfResultFile, croppedPriors, 0.1, OUT_DIR + fileNoExt + "/CausalPath-rank-based");
//			prepareCPDir(DATA_DIR + file, OUT_DIR + fileNoExt + "/CausalPath-classical");
		}

		for (String fileNoExt1 : rankedLists.keySet())
		{
			for (String fileNoExt2 : rankedLists.keySet())
			{
				String fileNoExt = fileNoExt1 + "_vs_" + fileNoExt2;

				FileUtil.mkdirs(OUT_DIR + "pairs/" + fileNoExt);

				List<String> preRankedList1 = rankedLists.get(fileNoExt1);
				List<String> preRankedList2 = rankedLists.get(fileNoExt2);

				List<String> rankedList1 = new ArrayList<>(preRankedList1);
				List<String> rankedList2 = new ArrayList<>(preRankedList2);
				rankedList1.retainAll(rankedList2);
				rankedList2.retainAll(rankedList1);

//				Map<String, Set<String>> geneSets = loadReactomeGeneSets();
//				reportDifferentialEnrichments(rankedList1, rankedList2, cropGeneSetsToSpecificList(geneSets, rankedList1), OUT_DIR + "pairs/" + fileNoExt + "/" + fileNoExt + ".reactome.xls");
//				geneSets = loadHallmarkGeneSets();
//				reportDifferentialEnrichments(rankedList1, rankedList2, cropGeneSetsToSpecificList(geneSets, rankedList1), OUT_DIR + "pairs/"  + fileNoExt + "/" +fileNoExt + ".hallmark.xls");

				preRankedList1 = signedRankedLists.get(fileNoExt1);
				preRankedList2 = signedRankedLists.get(fileNoExt2);
				rankedList1 = new ArrayList<>(preRankedList1);
				rankedList2 = new ArrayList<>(preRankedList2);
				rankedList1.retainAll(rankedList2);
				rankedList2.retainAll(rankedList1);

				Map<String, Map<String, Boolean>> croppedPriors = cropPriorsToSpecificList(priors, rankedList1);
				String tfResultFile = OUT_DIR + "pairs/" + fileNoExt + "/" + fileNoExt + ".TF-activity.xls";
				reportDifferentialTFTargetSignedEnricments(rankedList1, rankedList2, croppedPriors, tfResultFile);
			}
		}
	}

	private static List<String> readRankedList(String file)
	{
		return FileUtil.linesTabbedSkip(file, 1).map(t -> t[0]).collect(Collectors.toList());
	}

	private static List<String> readSignedRankedList(String file)
	{
		String[] header = FileUtil.readHeader(file);
		int fcInd = ArrayUtil.indexOf(header, "avg_log2FC");
		List<String> list = FileUtil.linesTabbedSkip(file, 1).filter(t -> !t[fcInd].startsWith("-")).map(t -> t[0]).collect(Collectors.toList());
		List<String> negList = FileUtil.linesTabbedSkip(file, 1).filter(t -> t[fcInd].startsWith("-")).map(t -> t[0]).collect(Collectors.toList());
		Collections.reverse(negList);
		list.addAll(negList);
		return list;
	}



	private static List<String> convertToHuman(List<String> rankedList)
	{
		List<String> humanRankedList = new ArrayList<>();

		Set<String> repeaters = new HashSet<>();
		Set<String> memory = new HashSet<>();
		for (String gene : rankedList)
		{
			Set<String> human = MGI.get().getCorrespondingHumanSymbols(gene);
			if (human.size() == 1)
			{
				String hGene = human.iterator().next();

				if (memory.contains(hGene))
				{
					System.err.println("Human gene repeats: " + gene + ", " + hGene);
					repeaters.add(hGene);
				}
				else
				{
					memory.add(hGene);
					humanRankedList.add(hGene);
				}
			}
		}

		humanRankedList.removeAll(repeaters);
		return humanRankedList;
	}

	private static Map<String, Set<String>> loadReactomeGeneSets()
	{
		return MSigDB.get().getSetsNameFiltered(name -> name.contains("REACTOME"));
	}

	private static Map<String, Set<String>> loadHallmarkGeneSets()
	{
		return MSigDB.get().getSetsNameFiltered(name -> name.contains("HALLMARK"));
	}

	private static void reportEnrichments(List<String> rankedList, Map<String, Set<String>> geneSets, String outFile) throws IOException
	{
		RankedListEnrichment.reportEnrichment(rankedList, geneSets, RAND_ITER, outFile);
	}

	private static void reportDifferentialEnrichments(List<String> rankedList1, List<String> rankedList2, Map<String, Set<String>> geneSets, String outFile) throws IOException
	{
		RankedListDifferentialEnrichment.reportEnrichment(rankedList1, rankedList2, geneSets, RAND_ITER, outFile);
	}

	private static void reportTFTargetSignedEnricments(List<String> rankedList, Map<String, Map<String, Boolean>> priors, String outFile) throws IOException
	{
		RankedListSignedEnrichment.reportEnrichment(rankedList, priors, RAND_ITER, outFile);
	}

	private static void reportDifferentialTFTargetSignedEnricments(List<String> rankedList1, List<String> rankedList2, Map<String, Map<String, Boolean>> priors, String outFile) throws IOException
	{
		RankedListSignedDifferentialEnrichment.reportEnrichment(rankedList1, rankedList2, priors, RAND_ITER, outFile);
	}

	private static Map<String, Set<String>> cropGeneSetsToSpecificList(Map<String, Set<String>> origSets, List<String> rankedList)
	{
		Set<String> set = new HashSet<>(rankedList);

		Map<String, Set<String>> geneSets = new HashMap<>();

		origSets.forEach((name, gset) ->
		{
			Set<String> croppedSet = new HashSet<>(gset);
			croppedSet.retainAll(set);
			if (croppedSet.size() > 3) geneSets.put(name, croppedSet);
		});

		return mergeIdenticalSets(geneSets);
	}

	private static Map<String, Map<String, Boolean>> cropPriorsToSpecificList(Map<String, Map<String, Boolean>> origPriors, List<String> rankedList)
	{
		Set<String> set = new HashSet<>(rankedList);

		Map<String, Map<String, Boolean>> priors = new HashMap<>();

		origPriors.forEach((name, gset) ->
		{
			Map<String, Boolean> croppedSet = new HashMap<>(gset);
			new HashSet<>(croppedSet.keySet()).stream().filter(gene -> !set.contains(gene)).forEach(croppedSet::remove);
			if (croppedSet.size() > 3) priors.put(name, croppedSet);
		});

		return mergeIdenticalMaps(priors);
	}

	private static Map<String, Set<String>> mergeIdenticalSets(Map<String, Set<String>> origSets)
	{
		Map<Set<String>, Set<String>> setToNames = new HashMap<>();

		origSets.forEach((name, set) ->
		{
			if (!setToNames.containsKey(set)) setToNames.put(set, new HashSet<>());
			setToNames.get(set).add(name);
		});

		Map<String, Set<String>> geneSets = new HashMap<>();
		setToNames.forEach((set, names) -> geneSets.put(CollectionUtil.merge(names, ", "), set));
		return geneSets;
	}

	private static Map<String, Map<String, Boolean>> mergeIdenticalMaps(Map<String, Map<String, Boolean>> origPriors)
	{
		Map<Map<String, Boolean>, Set<String>> mapToNames = new HashMap<>();

		origPriors.forEach((name, map) ->
		{
			if (!mapToNames.containsKey(map)) mapToNames.put(map, new HashSet<>());
			mapToNames.get(map).add(name);
		});

		Map<String, Map<String, Boolean>> priors = new HashMap<>();
		mapToNames.forEach((map, names) -> priors.put(CollectionUtil.merge(names, ", "), map));
		return priors;
	}

	private static Map<String, Map<String, Boolean>> readTFTargetPriors()
	{
		Map<String, Map<String, Boolean>> sets = new HashMap<>();

		FileUtil.linesTabbed("/home/ozgunbabur/Data/causal-priors.txt").filter(t -> t[1].endsWith("-expression")).forEach(t ->
		{
			if (!sets.containsKey(t[0])) sets.put(t[0], new HashMap<>());
			sets.get(t[0]).put(t[2], t[1].startsWith("up"));
		});

		return sets;
	}

	private static String getNodeColors(String file)
	{
		String[] header = FileUtil.readHeader(file);
		int idInd = 0;
		int fcInd = ArrayUtil.indexOf(header, "avg_log2FC");
		int pInd = ArrayUtil.indexOf(header, "p_val_adj");

		StringBuilder sb = new StringBuilder();
		ValToColor vtc = new ValToColor(new double[]{-10, 0, 10}, new Color[]{new Color(100, 100, 255), Color.WHITE, new Color(255, 100, 100)});

		Map<String, String> mColors = new HashMap<>();
		Map<String, String> mFC = new HashMap<>();
		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			String id = t[idInd];
			boolean neg = t[fcInd].startsWith("-");
			double p = Double.parseDouble(t[pInd]);
			if (p == 0) p = 1e-300;

			String color = vtc.getColorInString(-Math.log(p) * (neg ? -1 : 1));
			mColors.put(id, color);
			mFC.put(id, t[fcInd]);
		});

		Map<String, String> hColors = mColors.keySet().stream()
			.filter(mID -> MGI.get().getCorrespondingHumanSymbols(mID).size() == 1)
			.collect(Collectors.toMap(mID -> MGI.get().getCorrespondingHumanSymbols(mID).iterator().next(), mColors::get, (s1, s2) -> s1));

		Map<String, String> hFC = mFC.keySet().stream()
			.filter(mID -> MGI.get().getCorrespondingHumanSymbols(mID).size() == 1)
			.collect(Collectors.toMap(mID -> MGI.get().getCorrespondingHumanSymbols(mID).iterator().next(), mFC::get, (s1, s2) -> s1));

		sb.append("node\tall-nodes\tcolor\t255 255 255\nnode\tall-nodes\tbordercolor\t50 50 50");
		hColors.keySet().forEach(hID ->
			sb.append("\nnode\t").append(hID).append("\trppasite\t").append(hID).append("-rna|r|")
				.append(hColors.get(hID)).append("|50 50 50|").append(hFC.get(hID)));
		return sb.toString();
	}

	private static void prepareCPDir(String dataFile, String cpDir) throws IOException
	{
		String[] header = FileUtil.readHeader(dataFile);
		int idInd = 0;
		int fcInd = ArrayUtil.indexOf(header, "avg_log2FC");
		int pInd = ArrayUtil.indexOf(header, "p_val_adj");

		FileUtil.mkdirs(cpDir);
		BufferedWriter writer = FileUtil.newBufferedWriter(cpDir + "/data.tsv");
		writer.write("ID\tSymbols\tSites\tFeature\tEffect\tSignedP");

		FileUtil.linesTabbedSkip1(dataFile).forEach(t ->
		{
			String sym = t[idInd];
			Set<String> humanSymbols = MGI.get().getCorrespondingHumanSymbols(sym);
			if (humanSymbols.size() == 1)
			{
				sym = humanSymbols.iterator().next();
				boolean neg = t[fcInd].startsWith("-");
				double p = Double.parseDouble(t[pInd]);
				if (p == 0) p = 1e-300;
				if (neg) p = -p;
				FileUtil.lnwrite(sym + "-rna\t" + sym + "\t\tR\t\t" + p, writer);
			}
		});

		writer.close();

		FileUtil.writeStringToFile(PARAMETERS, cpDir + "/parameters.txt");
	}

	private static void generateNetworkResults(String dataFile, String tfResultFile, Map<String, Map<String, Boolean>> priors, double fdrThr, String outDir) throws IOException
	{
		String sifName = outDir + "/" + tfResultFile.substring(tfResultFile.lastIndexOf("/") + 1, tfResultFile.lastIndexOf("."));
		Optional<Double> optVal = FileUtil.linesTabbedSkip1(tfResultFile).filter(t -> Double.parseDouble(t[3]) <= fdrThr)
			.map(t -> Double.parseDouble(t[2])).max(Double::compareTo);

		if (optVal.isPresent())
		{
			FileUtil.mkdirs(outDir);
			Double pThr = optVal.get();

			BufferedWriter sifWriter = FileUtil.newBufferedWriter(sifName + ".sif");
			Map<String, Double> pvals = new HashMap<>();
			FileUtil.linesTabbedSkip1(tfResultFile).filter(t -> Double.parseDouble(t[2]) <= pThr).forEach(t ->
			{
				Set<String> targets = new HashSet<>(Arrays.asList(t[4].split(" ")));
				for (String tf : t[0].split(", "))
				{
					double p = Double.parseDouble(t[2]);
					if (p == 0) p = 1D / RAND_ITER;
					if (t[1].equals("-")) p = -p;
					pvals.put(tf, p);

					targets.forEach(target ->
						FileUtil.writeln(tf + "\t" + (priors.get(t[0]).get(target) ? "up" : "down") +
							"regulates-expression\t" + target, sifWriter));
				}
			});

			sifWriter.close();

			BufferedWriter fmtWriter = FileUtil.newBufferedWriter(sifName + ".format");
			FileUtil.write(getNodeColors(dataFile), fmtWriter);

			ValToColor vtc = new ValToColor(new double[]{Math.log(0.00001), 0, -Math.log(0.00001)}, new Color[]{new Color(100, 100, 255), Color.WHITE, new Color(255, 100, 100)});
			pvals.forEach((tf, p) ->
			{
				String letter = p > 0 ? "!" : "i";
				String boderColor = p > 0 ? "0 180 20" : "180 0 20";
				String featID = tf + (p > 0 ? "-activated" : "-inactivated");
				FileUtil.lnwrite("node\t" + tf + "\trppasite\t" + featID + "|" + letter + "|" + vtc.getColorInString(Math.log(Math.abs(p)) * (p > 0 ? -1 : 1)) + "|" + boderColor + "|" + Math.abs(p), fmtWriter);
				FileUtil.lnwrite("node\t" + tf + "\tbordercolor\t" + boderColor, fmtWriter);
			});

			fmtWriter.close();

			GraphWriter.convertSIFToJSON(sifName + ".sif", sifName + ".format", sifName.substring(0, sifName.lastIndexOf("/") + 1) + "causative.json");
		}
	}

	static String PARAMETERS = "proteomics-values-file = data.tsv\n" +
		"id-column = ID\n" +
		"symbols-column = Symbols\n" +
		"sites-column = Sites\n" +
		"feature-column = Feature\n" +
		"effect-column = Effect\n" +
		"value-transformation = signed-p-values\n" +
		"value-column = SignedP\n" +
		"threshold-for-data-significance = 0.1 rna\n" +
		"data-type-for-expressional-targets = rna\n" +
		"calculate-network-significance = true\n" +
		"permutations-for-significance = 10000\n" +
		"fdr-threshold-for-network-significance = 0.1\n" +
		"use-network-significance-for-causal-reasoning = true\n" +
		"minimum-potential-targets-to-consider-for-downstream-significance = 5\n" +
		"color-saturation-value = 10\n" +
		"show-all-genes-with-proteomic-data = false";
}
