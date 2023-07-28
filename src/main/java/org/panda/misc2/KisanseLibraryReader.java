package org.panda.misc2;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.Histogram;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class KisanseLibraryReader
{
	static final int PERCENT_THRESHOLD = 90;
	static final int RANK_THRESHOLD = 15;
	static final String DIR = "/home/ozgunbabur/Data/KinaseLibrary/";
	static final String SUPP_PATH = DIR + "KL_scores.csv";
	static final String SIF_PATH = DIR + "kinase-library-p"+ PERCENT_THRESHOLD +"-r" + RANK_THRESHOLD + ".sif";

	public static void main(String[] args) throws IOException
	{
//		load(PERCENT_THRESHOLD, RANK_THRESHOLD);
		convertToSIF(PERCENT_THRESHOLD, RANK_THRESHOLD);
//		compareCPAndKL();
	}

	private static Map<String, Map<String, Map<String, double[]>>> load(double percentileThr, int rankThr)
	{
		Map<String, Map<String, Map<String, double[]>>> map = new HashMap<>();

		String[] header = FileUtil.readHeader(SUPP_PATH);
		int geneInd = ArrayUtil.indexOf(header, "Gene");
		int siteInd = ArrayUtil.indexOf(header, "Phosphosite");
		int kinStartInd = ArrayUtil.indexOf(header, "AAK1_percentile");

		int[] kinCnt = new int[212];
		Map<String, Integer> subsCnt = new HashMap<>();

		FileUtil.linesTabbedSkip1(SUPP_PATH)
			.filter(t -> !t[geneInd].isEmpty() && !t[geneInd].contains(";"))
			.forEach(t ->
		{
			String target = t[geneInd];
			String site = t[siteInd];

			int cnt = 0;

			for (int i = kinStartInd; i < t.length; i+=2)
			{
				double percentile = Double.parseDouble(t[i]);

				if (percentile >= percentileThr)
				{
					cnt++;

					String kinase = header[i].substring(0, header[i].lastIndexOf("_"));
					subsCnt.put(kinase, subsCnt.getOrDefault(kinase, 0) + 1);
				}
			}

			kinCnt[cnt]++;

//			if (cnt <= rankThr)
			for (int i = kinStartInd; i < t.length; i += 2)
			{
				String kinase = header[i].substring(0, header[i].lastIndexOf("_"));
				double percentile = Double.parseDouble(t[i]);
				int rank = Integer.parseInt(t[i+1]);

				if (percentile >= percentileThr && rank <= rankThr)
				{
					cnt++;
					if (!map.containsKey(kinase)) map.put(kinase, new HashMap<>());

					Map<String, Map<String, double[]>> targetMap = map.get(kinase);
					if (!targetMap.containsKey(target)) targetMap.put(target, new HashMap<>());
					if (!targetMap.get(target).containsKey(site)) targetMap.get(target).put(site, new double[]{percentile, rank});
				}
			}
		});

		double total = ArrayUtil.sum(kinCnt);
		int[] cum = new int[kinCnt.length];
		cum[0] = kinCnt[0];
		for (int i = 1; i < kinCnt.length; i++)
		{
			cum[i] = cum[i-1] + kinCnt[i];
		}
		for (int i = 0; i < cum.length; i++)
		{
//			System.out.println(i + "\t" + (cum[i] / total));
//			System.out.println(i + "\t" + kinCnt[i]);
		}

//		Histogram h = new Histogram(100);
//		h.setBorderAtZero(true);
//		subsCnt.forEach((s, integer) -> h.count(integer));
//		h.print();

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
