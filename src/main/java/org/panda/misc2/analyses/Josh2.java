package org.panda.misc2.analyses;

import org.panda.misc2.TFEnrichment;
import org.panda.resource.MGI;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.RankedListSignedEnrichment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class Josh2
{
//	public static final String DATA_DIR = "/home/ozgunbabur/Data/Josh/take3/DEGs/";
	public static final String DATA_DIR = "/home/ozgunbabur/Data/Josh/take3/Shailja-alternative/";
//	public static final String ANALYSIS_DIR = "/home/ozgunbabur/Analyses/Josh/take3-CP-redo/";
	public static final String ANALYSIS_DIR = "/home/ozgunbabur/Analyses/Josh/Shailja-alternative/";

	public static void main(String[] args) throws IOException
	{
//		batchConvertMouseToHuman(DATA_DIR, ANALYSIS_DIR);
//		runRankBasedEnrichment(ANALYSIS_DIR);
		addRankBasedEnrichmentToCausalPathRecursive(ANALYSIS_DIR, 0.1);
//		temp();
	}


	static Map<String, Integer> loadClusters()
	{
		Map<String, String> ss = FileUtil.readMap(DATA_DIR + "meta_data.tsv", "\t", "cell_ID", "seurat_clusters");
		return ss.keySet().stream().collect(Collectors.toMap(c -> c, c -> Integer.parseInt(ss.get(c))));
	}

	static Map<String, double[]> readCellCoordinates()
	{
		Map<String, String> xMap = FileUtil.readMap(DATA_DIR + "umap_coords.tsv", "\t", "cell_ID", "UMAP_1");
		Map<String, String> yMap = FileUtil.readMap(DATA_DIR + "umap_coords.tsv", "\t", "cell_ID", "UMAP_2");
		return xMap.keySet().stream().collect(Collectors.toMap(c -> c, c -> new double[]{Double.parseDouble(xMap.get(c)), Double.parseDouble(yMap.get(c))}));
	}

	static void printLocDist(int cluster)
	{
		Map<String, Integer> clusMap = loadClusters();
		Map<String, double[]> coorMap = readCellCoordinates();

		double res = 0.1;
		Histogram hx = new Histogram(res);
		Histogram hy = new Histogram(res);
		clusMap.forEach((cell, clus) ->
		{
			if (clus == cluster)
			{
				double[] coor = coorMap.get(cell);
				hx.count(coor[0]);
				hy.count(coor[1]);
			}
		});
		hx.print();
		hy.print();
	}

	static void printCellsInRange(int clust, double xMin, double xMax, double yMin, double yMax)
	{
		Map<String, Integer> clusMap = loadClusters();
		Map<String, double[]> coorMap = readCellCoordinates();

		clusMap.forEach((cell, clus) ->
		{
			if (clust >= 0 && clus == clust)
			{
				double[] crd = coorMap.get(cell);

				if (crd[0] >= xMin && crd[0] <= xMax && crd[1] >= yMin && crd[1] <= yMax)
				{
					System.out.println(cell + "\t" + crd[0] + "\t" + crd[1]);
				}
			}
			else
			{
				double[] crd = coorMap.get(cell);

				if (crd[0] >= xMin && crd[0] <= xMax && crd[1] >= yMin && crd[1] <= yMax)
				{
					System.out.println(cell + "\t" + clus + "\t" + crd[0] + "\t" + crd[1]);
				}
			}
		});
	}

	static Map<Integer, String> getClusterRepresentatives()
	{
		Map<String, Integer> clusMap = loadClusters();
		Map<String, double[]> coorMap = readCellCoordinates();

		Map<Integer, ArrayList> xLists = new HashMap<>();
		Map<Integer, ArrayList> yLists = new HashMap<>();
		for (int i = 0; i < 12; i++)
		{
			xLists.put(i, new ArrayList());
			yLists.put(i, new ArrayList());
		}

		clusMap.forEach((cell, clus) ->
		{
			double[] crd = coorMap.get(cell);
			xLists.get(clus).add(crd[0]);
			yLists.get(clus).add(crd[1]);
		});

		Map<Integer, double[]> centers = new HashMap<>();
		xLists.forEach((clus, xList) ->
		{
			ArrayList yList = yLists.get(clus);
			double[] coor = new double[]{CollectionUtil.averageInListOfDouble(xList), CollectionUtil.averageInListOfDouble(yList)};
			centers.put(clus, coor);
		});

		Map<Integer, String> repMap = new HashMap<>();
		Map<Integer, Double> disMap = new HashMap<>();

		clusMap.forEach((cell, clus) ->
		{
			double[] loc = coorMap.get(cell);
			double[] ctr = centers.get(clus);
			double dist = ((loc[0] - ctr[0]) * (loc[0] - ctr[0])) + ((loc[1] - ctr[1]) * (loc[1] - ctr[1]));

			if (!repMap.containsKey(clus) || disMap.get(clus) > dist)
			{
				repMap.put(clus, cell);
				disMap.put(clus, dist);
			}
		});

		repMap.forEach((cl, cell) -> System.out.println(cl + "\t" + cell + "\t" + Arrays.toString(coorMap.get(cell))));
		return repMap;
	}

	static void logTransformMatrix() throws IOException
	{
		String inFile = DATA_DIR + "expression_matrix.tsv";
		String firstLine = FileUtil.lines(inFile).findFirst().get();

		BufferedWriter writer = FileUtil.newBufferedWriter(DATA_DIR + "expression_matrix_log_transformed.tsv");
		writer.write(firstLine);

		FileUtil.linesTabbedSkip1(inFile).forEach(t ->
		{
			String gene = t[0];
			FileUtil.lnwrite(gene, writer);
			for (int i = 1; i < t.length; i++)
			{
				double v = Double.parseDouble(t[i]);
				v = Math.log1p(v);
				FileUtil.tab_write(v, writer);
			}
		});

		writer.close();
	}

	static void batchConvertMouseToHuman(String inDir, String outDir) throws IOException
	{
		for (File file : new File(inDir).listFiles())
		{
			if (file.getName().endsWith(".csv"))
			{
				String name = file.getName();
				name = name.substring(0, name.lastIndexOf("."));

				String outFile = outDir + name + "/data.tsv";
				FileUtil.mkdirsOfFilePath(outFile);
				convertMouseToHuman(file.getPath(), outFile);
			}
		}
	}

	static void convertMouseToHuman(String inFile, String outFile) throws IOException
	{
		String[] header = FileUtil.readHeader(inFile);
		int pInd = ArrayUtil.indexOf(header, "p_val");
		int fcInd = ArrayUtil.indexOf(header, "avg_log2FC");

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("ID\tSymbols\tSites\tFeature\tEffect\tSignedP");
		Map<String, String> lineMap = new HashMap<>();
		Set<String> repeaters = new HashSet<>();
		Set<String> memory = new HashSet<>();
		FileUtil.linesTabbedSkip1(inFile).forEach(t ->
		{
			if (t.length < 2 || t[1].isEmpty()) return;

			String mSym = t[0];

			Set<String> human = MGI.get().getCorrespondingHumanSymbols(mSym);
			if (human.size() == 1)
			{
				String hGene = human.iterator().next();

				if (memory.contains(hGene))
				{
					System.err.println("Human gene repeats: " + mSym + ", " + hGene);
					repeaters.add(hGene);
				}
				else
				{
					memory.add(hGene);
					double p = Double.parseDouble(t[pInd]);
					if (p == 0)
					{
						p = 1e-323;
					}
					if (t[fcInd].startsWith("-"))
					{
						p = -p;
					}
					String line = hGene + "-rna\t" + hGene + "\t\tR\t\t" + p;
					lineMap.put(hGene, line);
				}
			}
		});
		repeaters.forEach(lineMap::remove);
		lineMap.values().forEach(l -> FileUtil.lnwrite(l, writer));
		writer.close();
	}

	static void temp()
	{
		Set<String> targets = FileUtil.linesTabbed("/home/ozgunbabur/Documents/causal-priors.txt").filter(t -> t[0].equals("IRF1") &&
			(t[1].equals("upregulates-expression") || t[1].equals("downregulates-expression"))).map(t -> t[2]).collect(Collectors.toSet());

		System.out.println("targets = " + targets.size());

		Set<String> oldSet = FileUtil.getTermsInTabDelimitedColumn("/home/ozgunbabur/Data/Josh/take3/DEGs/RS_LPvsLP.xls", 0, 1);
		Set<String> newSet = FileUtil.getTermsInTabDelimitedColumn("/home/ozgunbabur/Data/Josh/take3/BPK2022vS4/RS_LPvLP.xls", 0, 1);

//		Set<String> oldSet = FileUtil.getTermsInTabDelimitedColumn("/home/ozgunbabur/Analyses/Josh/take3-CP-redo/RS_LPvsLP/data.tsv", 1, 1);
//		Set<String> newSet = FileUtil.getTermsInTabDelimitedColumn("/home/ozgunbabur/Analyses/Josh/BPK2022vS4/RS_LPvLP/data.tsv", 1, 1);
//		oldSet.retainAll(targets);
//		newSet.retainAll(targets);

		CollectionUtil.printVennCounts(oldSet, newSet);

		System.out.println(" ========== ");

		Map<String, Boolean> map1 = FileUtil.linesTabbedSkip1("/home/ozgunbabur/Data/Josh/take3/DEGs/RS_LPvsLP.xls").filter(t ->      Double.parseDouble(t[5]) < 0.000000000000001).collect(Collectors.toMap(t -> t[0], t -> !t[2].startsWith("-")));
		Map<String, Boolean> map2 = FileUtil.linesTabbedSkip1("/home/ozgunbabur/Data/Josh/take3/BPK2022vS4/RS_LPvLP.xls").filter(t -> Double.parseDouble(t[5]) < 0.000000000000001).collect(Collectors.toMap(t -> t[0], t -> !t[2].startsWith("-")));

		int[] cnt = new int[]{0, 0};
		map1.keySet().stream().filter(map2.keySet()::contains).forEach(g ->
		{
			if (map1.get(g).equals(map2.get(g))) cnt[0]++;
			else
			{
				cnt[1]++;
//				System.out.println(g);
			}
		});

		System.out.println("cnt = " + Arrays.toString(cnt));

	}

	static void batchConvertToHuman() throws IOException
	{
		String outBase = "/home/ozgunbabur/Analyses/Josh/groups/";
		String inBase = "/home/ozgunbabur/Data/Josh/diffexp/groups/";

		for (File file : new File(inBase).listFiles())
		{
			String name = file.getName();
			name = name.substring(0, name.lastIndexOf("."));
			FileUtil.mkdirs(outBase + name);
			convertMouseToHuman(file.getPath(), outBase + name + "/data.tsv");
			FileUtil.copyFile(inBase + "../../parameters.txt", outBase + name + "/parameters.txt") ;
		}
	}

	static void runRankBasedEnrichment(String base) throws IOException
	{
		for (File dir : new File(base).listFiles())
		{
			System.out.println("dir = " + dir);
			Map<String, Map<String, Boolean>> rawPriors = TFEnrichment.readPriors();
			List<String> ranking = TFEnrichment.readRankedIDsFromCPFile(dir + "/data.tsv", "R");
			Map<String, Map<String, Boolean>> priors = TFEnrichment.convertPriors(rawPriors, ranking, 3);
			RankedListSignedEnrichment.reportEnrichment(ranking, priors, 1000000, dir + "/tf-enrichment.tsv");
		}
	}

	static void addRankBasedEnrichmentToCausalPathRecursive(String base, double fdrThr) throws IOException
	{
		FileUtil.processDirsRecursive(new File(base), dir ->
		{
			if (FileUtil.exists(dir + "/parameters.txt") && FileUtil.exists(dir + "/tf-enrichment.tsv"))
			{
				addRankBasedEnrichmentToCausalPath(dir.getPath(), fdrThr);
			}
		});
	}

	static void addRankBasedEnrichmentToCausalPath(String dir, double fdrThr)
	{
		List<String[]> tokenList = FileUtil.linesTabbedSkip1(dir + "/tf-enrichment.tsv").collect(Collectors.toList());
		int k = -1;
		for (int i = 0; i < tokenList.size(); i++)
		{
			String[] t = tokenList.get(i);
			if (Double.parseDouble(t[3]) < fdrThr) k = i;
		}

		if (k < 0) return;

		StringBuilder sb = new StringBuilder();
		for (int i = 0; i <= k; i++)
		{
			String[] t = tokenList.get(i);
			sb.append("\ngene-activity = ").append(t[0]).append(" ").append(t[1]);
		}
		FileUtil.addStringToFile(sb.toString(), dir + "/parameters.txt");
//		System.out.println("sb.toString() = " + sb.toString());
	}
}
