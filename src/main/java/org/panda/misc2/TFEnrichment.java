package org.panda.misc2;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.RankedListSignedEnrichment;
import org.panda.utility.statistics.RankedListSignedGroupedEnrichment;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class TFEnrichment
{
	public static void main(String[] args) throws IOException
	{
		Map<String, Map<String, Boolean>> rawPriors = readPriors();

		String dataFile = "/home/ozgunbabur/Data/Aslan/gpvi_mass_spec/C2S.tsv";
		List<String> rankedList = readRankedIDsFromCPFile(dataFile, "R");

		// for sanity check
//		Collections.shuffle(rankedList);

		Map<String, Map<String, Boolean>> priors = convertPriors(rawPriors, rankedList, 3);

		RankedListSignedEnrichment.reportEnrichment(rankedList, priors, 1000000, dataFile.substring(0, dataFile.lastIndexOf(".")) + "-tf-enrichment.tsv");
	}

	public static List<String> readRankedIDsFromCPFile(String file, String feature)
	{
		String[] header = FileUtil.readHeader(file);
		int idInd = ArrayUtil.indexOf(header, "Symbols");
		int pInd = ArrayUtil.indexOf(header, "SignedP");
		int featInd = ArrayUtil.indexOf(header, "Feature", "Modification");
		Map<String, Double> idToP = FileUtil.linesTabbedSkip1(file).filter(t -> t[featInd].equals(feature)).collect(Collectors.toMap(t -> t[idInd], t -> Double.parseDouble(t[pInd]), (d1, d2) -> d1));
		List<String> ids = new ArrayList<>(idToP.keySet());
		ids.sort((id1, id2) ->
		{
			double p1 = idToP.get(id1);
			double p2 = idToP.get(id2);
			if (Math.signum(p1) == Math.signum(p2)) return Double.compare(p1, p2);
			else return Double.compare(p2, p1);
		});
		return ids;
	}

	public static Map<String, Map<String, Boolean>> convertPriors(Map<String, Map<String, Boolean>> priors, List<String> consider, int minTarg)
	{
		Set<String> keep = new HashSet<>(consider);
		Map<String, Map<String, Boolean>> converted = new HashMap<>();

		for (String tf : priors.keySet())
		{
			Map<String, Boolean> targets = new HashMap<>();
			converted.put(tf, targets);

			for (String target : priors.get(tf).keySet())
			{
				if (keep.contains(target))
				{
					boolean sign = priors.get(tf).get(target);
					targets.put(target, sign);
				}
			}
		}

		Set<String> remove = new HashSet<>();
		converted.keySet().stream().filter(tf -> converted.get(tf).size() < minTarg).forEach(remove::add);
		remove.forEach(converted::remove);

		return mergeIdenticalSets(converted);
	}

	private static Map<String, Map<String, Boolean>> mergeIdenticalSets(Map<String, Map<String, Boolean>> prePriors)
	{
		Map<Map<String, Boolean>, Set<String>> mapsToNames = new HashMap<>();
		prePriors.forEach((name, maps) ->
		{
			if (!mapsToNames.containsKey(maps)) mapsToNames.put(maps, new HashSet<>());
			mapsToNames.get(maps).add(name);
		});

		Map<String, Map<String, Boolean>> priors = new HashMap<>();
		mapsToNames.forEach((maps, names) -> priors.put(CollectionUtil.merge(names, ", "), maps));
		return priors;
	}

	public static Map<String, Map<String, Boolean>> readPriors()
	{
		Map<String, Map<String, Boolean>> priors = new HashMap<>();

		FileUtil.linesTabbed("/home/ozgunbabur/Data/causal-priors.txt")
			.filter(t -> t.length > 2 && t[1].contains("expression")).forEach(t ->
		{
			if (!priors.containsKey(t[0])) priors.put(t[0], new HashMap<>());
			Map<String, Boolean> targets = priors.get(t[0]);
			boolean direction = t[1].startsWith("u");
			targets.put(t[2], direction);
		});
		return priors;
	}
}
