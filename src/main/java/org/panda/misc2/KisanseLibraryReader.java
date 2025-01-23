package org.panda.misc2;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class KisanseLibraryReader
{
	static final int PERCENT_THRESHOLD = 90;
//	static final int RANK_THRESHOLD = 15;
	static final int RANK_THRESHOLD = 5;
	static final String DIR = "/home/ozgunbabur/Data/KinaseLibrary/";
	static final String SCORES_ST_PATH = DIR + "KL_scores_ST.csv";
	static final String SCORES_Y_PATH = DIR + "KL_scores_Y.csv";
	static final String SIF_PATH = DIR + "kinase-library-p"+ PERCENT_THRESHOLD +"-r" + RANK_THRESHOLD + ".sif";

	public static void main(String[] args) throws IOException
	{
//		load(PERCENT_THRESHOLD, RANK_THRESHOLD);
		convertToSIF(PERCENT_THRESHOLD, RANK_THRESHOLD);
//		compareCPAndKL();
//		compareKinases();
	}

	public static void compareKinases()
	{
		String[] t1 = FileUtil.readFirstLine(SCORES_ST_PATH).split("\t");
		String[] t2 = FileUtil.readFirstLine(SCORES_Y_PATH).split("\t");

		Set<String> s1 = new HashSet<>();
		Set<String> s2 = new HashSet<>();

		for (String s : t1)
		{
			if (s.endsWith("_rank")) s1.add(s.substring(0, s.length() - 5));
		}
		for (String s : t2)
		{
			if (s.endsWith("_rank")) s2.add(s.substring(0, s.length() - 5));
		}
		CollectionUtil.printVennSets(s1, s2);
	}

	private static Map<String, Map<String, Map<String, double[]>>> load(double percentileThr, int rankThr)
	{
		Map<String, Map<String, Map<String, double[]>>> map = new HashMap<>();

		String[] headerST = FileUtil.readHeader(SCORES_ST_PATH);
		int geneIndST = ArrayUtil.indexOf(headerST, "Gene");
		int siteIndST = ArrayUtil.indexOf(headerST, "Phosphosite");
		int kinStartIndST = ArrayUtil.indexOf(headerST, "AAK1_percentile");

		Map<String, Integer> subsCnt = new HashMap<>();

		FileUtil.linesTabbedSkip1(SCORES_ST_PATH)
			.filter(t -> !t[geneIndST].isEmpty() && !t[geneIndST].contains(";"))
			.forEach(t ->
		{
			String target = t[geneIndST];
			String site = t[siteIndST];

			for (int i = kinStartIndST; i < t.length; i+=2)
			{
				double percentile = Double.parseDouble(t[i]);

				if (percentile >= percentileThr)
				{
					String kinase = headerST[i].substring(0, headerST[i].lastIndexOf("_"));
					subsCnt.put(kinase, subsCnt.getOrDefault(kinase, 0) + 1);
				}
			}

			for (int i = kinStartIndST; i < t.length; i += 2)
			{
				String kinase = headerST[i].substring(0, headerST[i].lastIndexOf("_"));
				double percentile = Double.parseDouble(t[i]);
				int rank = Integer.parseInt(t[i+1]);

				if (percentile >= percentileThr && rank <= rankThr)
				{
					if (!map.containsKey(kinase)) map.put(kinase, new HashMap<>());

					Map<String, Map<String, double[]>> targetMap = map.get(kinase);
					if (!targetMap.containsKey(target)) targetMap.put(target, new HashMap<>());
					if (!targetMap.get(target).containsKey(site)) targetMap.get(target).put(site, new double[]{percentile, rank});
				}
			}
		});

		String[] headerY = FileUtil.readHeader(SCORES_Y_PATH);
		int geneIndY = ArrayUtil.indexOf(headerY, "Gene");
		int siteIndY = ArrayUtil.indexOf(headerY, "Phosphosite");
		int kinStartIndY = ArrayUtil.indexOf(headerY, "ABL_percentile");

		FileUtil.linesTabbedSkip1(SCORES_Y_PATH)
			.filter(t -> !t[geneIndY].isEmpty() && !t[geneIndY].contains(";"))
			.forEach(t ->
		{
			String target = t[geneIndY];
			String site = t[siteIndY];

			for (int i = kinStartIndY; i < t.length; i+=2)
			{
				double percentile = Double.parseDouble(t[i]);

				if (percentile >= percentileThr)
				{
					String kinase = headerY[i].substring(0, headerY[i].lastIndexOf("_"));
					subsCnt.put(kinase, subsCnt.getOrDefault(kinase, 0) + 1);
				}
			}

			for (int i = kinStartIndY; i < t.length; i += 2)
			{
				String kinase = headerY[i].substring(0, headerY[i].lastIndexOf("_"));
				double percentile = Double.parseDouble(t[i]);
				int rank = Integer.parseInt(t[i+1]);

				if (percentile >= percentileThr && rank <= rankThr)
				{
					if (!map.containsKey(kinase)) map.put(kinase, new HashMap<>());

					Map<String, Map<String, double[]>> targetMap = map.get(kinase);
					if (!targetMap.containsKey(target)) targetMap.put(target, new HashMap<>());
					if (!targetMap.get(target).containsKey(site)) targetMap.get(target).put(site, new double[]{percentile, rank});
				}
			}
		});

		return map;
	}

	private static void convertToSIF(double percentileThr, int upstKinaseThr) throws IOException
	{
		Map<String, Map<String, Map<String, double[]>>> map = load(percentileThr, upstKinaseThr);
		BufferedWriter writer = FileUtil.newBufferedWriter(SIF_PATH);

		List<String> kinases = map.keySet().stream().sorted().collect(Collectors.toList());

		for (String kinase : kinases)
		{
			Map<String, Map<String, double[]>> tMap = map.get(kinase);

			List<String> targets = tMap.keySet().stream().sorted().collect(Collectors.toList());

			for (String target : targets)
			{
				Map<String, double[]> siteMap = tMap.get(target);

				List<String> sites = new ArrayList<>(siteMap.keySet());
				sites.sort(Comparator.comparingInt(s -> Integer.parseInt(s.substring(1))));

				StringBuilder sb = new StringBuilder();
				boolean start = true;
				for (String site : sites)
				{
					double[] pr = siteMap.get(site);

					if (start) start = false;
					else sb.append("|");

					sb.append(pr[0]).append(";").append((int) pr[1]);
				}

				FileUtil.writeln(kinase + "\tphosphorylates\t" + target + "\t\t" + CollectionUtil.merge(sites, ";") + "\t" + sb, writer);
			};
		};

		writer.close();
	}

	private static void compareCPAndKL()
	{
		Set<String> kl = FileUtil.linesTabbed(SIF_PATH).filter(t -> t[1].equals("phosphorylates")).map(t -> t[0] + " " + t[2]).collect(Collectors.toSet());
		Set<String> cp = FileUtil.linesTabbed("/home/ozgunbabur/Documents/causal-priors.txt").filter(t -> t[1].equals("phosphorylates")).map(t -> t[0] + " " + t[2]).collect(Collectors.toSet());

		CollectionUtil.printVennCounts(kl, cp);
	}

	public static Map[] getScoresFromSpecialSIF(String sifFile)
	{
		Map<String, Double> percentageMap = new HashMap<>();
		Map<String, Integer> rankMap = new HashMap<>();

		FileUtil.linesTabbed(sifFile).forEach(t ->
		{
			String keyPrefix = t[0] + " " + t[2] + " ";

			String[] scores = t[5].split("\\|");
			String[] sites = t[4].split(";");
			for (int i = 0; i < sites.length; i++)
			{
				String key = keyPrefix + sites[i];
				String[] s = scores[i].split(";");
				percentageMap.put(key, Double.parseDouble(s[0]));
				rankMap.put(key, Integer.parseInt(s[1]));
			}
		});
	
		return new Map[]{percentageMap, rankMap};
	}
}
