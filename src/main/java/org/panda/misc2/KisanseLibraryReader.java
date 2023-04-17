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
	static final int UPSTREAM_KINASE_THRESHOLD = 20;
	static final String DIR = "/home/ozgunbabur/Data/KinaseLibrary/";
	static final String SUPP_PATH = DIR + "KL_scores.csv";
	static final String SIF_PATH = DIR + "kinase-library-"+ PERCENT_THRESHOLD +"-" + UPSTREAM_KINASE_THRESHOLD + ".sif";

	public static void main(String[] args) throws IOException
	{
//		load(90);
		convertToSIF(PERCENT_THRESHOLD, UPSTREAM_KINASE_THRESHOLD);
//		compareCPAndKL();
	}

	private static Map<String, Map<String, List<String>>> load(double percentileThr, int upstrKinaseThr)
	{
		Map<String, Map<String, List<String>>> map = new HashMap<>();

		String[] header = FileUtil.readHeader(SUPP_PATH);
		int geneInd = ArrayUtil.indexOf(header, "Gene");
		int siteInd = ArrayUtil.indexOf(header, "Phosphosite");
		int kinStartInd = ArrayUtil.indexOf(header, "AAK1_percentile");

		int[] kinCnt = new int[212];

		FileUtil.linesTabbedSkip1(SUPP_PATH).filter(t -> !t[geneInd].isEmpty() && !t[geneInd].contains(";")).forEach(t ->
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
				}
			}

			kinCnt[cnt]++;

			if (cnt <= upstrKinaseThr)
			{
				for (int i = kinStartInd; i < t.length; i += 2)
				{
					String kinase = header[i].substring(0, header[i].lastIndexOf("_"));
					double percentile = Double.parseDouble(t[i]);

					if (percentile >= percentileThr)
					{
						cnt++;
						if (!map.containsKey(kinase)) map.put(kinase, new HashMap<>());

						Map<String, List<String>> targetMap = map.get(kinase);
						if (!targetMap.containsKey(target)) targetMap.put(target, new ArrayList<>());
						if (!targetMap.get(target).contains(site)) targetMap.get(target).add(site);
					}
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
			System.out.println(i + "\t" + kinCnt[i]);
		}

		return map;
	}

	private static void convertToSIF(double percentileThr, int upstKinaseThr) throws IOException
	{
		Map<String, Map<String, List<String>>> map = load(percentileThr, upstKinaseThr);
		BufferedWriter writer = FileUtil.newBufferedWriter(SIF_PATH);

		map.forEach((kinase, tMap) ->
		{
			tMap.forEach((target, sites) ->
			{
				sites.sort(Comparator.comparingInt(s -> Integer.parseInt(s.substring(1))));
				FileUtil.writeln(kinase + "\tphosphorylates\t" + target + "\t\t" + CollectionUtil.merge(sites, ";"), writer);
			});
		});

		writer.close();
	}

	private static void compareCPAndKL()
	{
		Set<String> kl = FileUtil.linesTabbed(SIF_PATH).filter(t -> t[1].equals("phosphorylates")).map(t -> t[0] + " " + t[2]).collect(Collectors.toSet());
		Set<String> cp = FileUtil.linesTabbed("/home/ozgunbabur/Documents/causal-priors.txt").filter(t -> t[1].equals("phosphorylates")).map(t -> t[0] + " " + t[2]).collect(Collectors.toSet());

		CollectionUtil.printVennCounts(kl, cp);
	}
}
