package org.panda.misc2.analyses;

import org.panda.resource.MGIVertebrateHomology;
import org.panda.resource.MSigDB;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.RankedListEnrichmentAnalytical;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class Linden
{
	public static final String BASE = "/home/ozgunbabur/Data/Linden/LiverRNA/RNA seq/DESeq2/";
	public static final String EN_OUT = BASE + "/MSigDBEnrichments/";
	public static final String RANK_OUT = BASE + "/ranked-lists/";

	public static void main(String[] args) throws IOException
	{
		calcMSigDBEnrichments();
	}

	public static void calcMSigDBEnrichments() throws IOException
	{
		Map<String, Set<String>> geneSets = MSigDB.get().getGeneSets();
		for (File dir : new File(BASE).listFiles())
		{
			String file = dir.getPath() + "/DESeq2_output/DESeq2_stats_allGenes_CookTested.txt";
			if (FileUtil.exists(file))
			{
				List<String> rankedList = getTheRankedList(file);
				RankedListEnrichmentAnalytical.reportEnrichment(rankedList, geneSets, EN_OUT + dir.getName() + ".tsv");

				writeRatRank(file, dir.getName());
			}
		}
	}

	public static void writeRatRank(String file, String caseName)
	{
		Map<String, Double> rMap = FileUtil.linesTabbedSkip1(file).filter(t -> t[0].length() > 1 && !t[0].equals("NA"))
			.collect(Collectors.toMap(t -> t[0].replaceAll("\"", ""), t -> Double.parseDouble(t[6]), Math::min));

		BufferedWriter writer = FileUtil.newBufferedWriter(RANK_OUT + caseName + ".txt");
		rMap.keySet().stream().sorted(Comparator.comparing(rMap::get)).forEach(g -> FileUtil.writeln(g, writer));
		FileUtil.closeWriter(writer);
	}

	public static List<String> getTheRankedList(String file)
	{
		Map<String, Double> rMap = FileUtil.linesTabbedSkip1(file).filter(t -> t[0].length() > 2)
			.collect(Collectors.toMap(t -> t[0].replaceAll("\"", ""), t -> Double.parseDouble(t[6]), Math::min));

		Set<String> multiple = new HashSet<>();
		Map<String, Double> hMap = new HashMap<>();
		rMap.forEach((rSym, val) ->
		{
			Set<String> syms = MGIVertebrateHomology.get().getCorrespondingHumanSymbols(rSym, MGIVertebrateHomology.Organism.RAT);
			if (syms.size() == 1)
			{
				String hSym = syms.iterator().next();
				if (!hMap.containsKey(hSym))
				{
					hMap.put(hSym, val);
				}
				else
				{
					multiple.add(hSym);
				}
			}
			else
			{
				String hSym = rSym.toUpperCase();
				if (syms.contains(hSym))
				{
					if (!hMap.containsKey(hSym))
					{
						hMap.put(hSym, val);
					}
					else
					{
						multiple.add(hSym);
					}
				}
			}
		});

		multiple.forEach(hMap::remove);

		List<String> list = new ArrayList<>(hMap.keySet());
		list.sort(Comparator.comparing(hMap::get));
		return list;
	}
}
