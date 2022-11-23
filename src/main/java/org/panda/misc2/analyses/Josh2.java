package org.panda.misc2.analyses;

import org.panda.misc2.TFEnrichment;
import org.panda.resource.MGI;
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
	public static final String DATA_DIR = "/home/ozgunbabur/Data/Josh/";
	public static final String ANALYSIS_DIR = "/home/ozgunbabur/Analyses/Josh/trajectories/";

	public static void main(String[] args) throws IOException
	{
//		printLocDist(6);
//		printCellsInRange(-1, -5.1, -4.9, 4.9, 5.1);
//		getClusterRepresentatives();
//		logTransformMatrix();
//		convertMouseToHuman(DATA_DIR + "diffexp.tsv", ANALYSIS_DIR + "data.tsv");
//		batchConvertToHuman();
//		runRankBasedEnrichment();
//		addRankBasedEnrichmentToCausalPathRecursive(ANALYSIS_DIR, 0.1);
		temp();
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

	static void convertMouseToHuman(String inFile, String outFile) throws IOException
	{
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
					if (Double.parseDouble(t[1]) == 0)
					{
						Double p = 1e-323;
						if (t[1].startsWith("-")) p = -p;
						t[1] = String.valueOf(p);
					}
					String line = hGene + "-rna\t" + hGene + "\t\tR\t\t" + t[1];
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
		int[] max = new int[]{0};
		FileUtil.linesTabbedSkip1("/home/ozgunbabur/Data/Josh/diffexp/luminal-vs-basal.tsv").filter(t -> t.length > 1).forEach(t ->
		{
			if (t[1].substring(1).contains("-"))
			{
				int x = Integer.parseInt(t[1].substring(t[1].lastIndexOf("-") + 1));
				if (x > max[0]) max[0] = x;
			}
		});
		System.out.println("max[0] = " + max[0]);

		String s = "1e-323";
		double d = Double.parseDouble(s);
		System.out.println("d = " + d);
		System.out.println("String.valueOf(d) = " + String.valueOf(d));
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
			FileUtil.copyFile(inBase + "../parameters.txt", outBase + name + "/parameters.txt") ;
		}
	}

	static void runRankBasedEnrichment() throws IOException
	{
		String base = "/home/ozgunbabur/Analyses/Josh/trajectories/";

		for (File dir : new File(base).listFiles())
		{
			System.out.println("dir = " + dir);
			Map<String, Map<String, Boolean>> rawPriors = TFEnrichment.readPriors();
			List<String> ranking = TFEnrichment.readRankedIDsFromCPFile(dir + "/data.tsv", "R");
			Map<String, Map<String, Boolean>> priors = TFEnrichment.convertPriors(rawPriors, ranking, 5);
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
