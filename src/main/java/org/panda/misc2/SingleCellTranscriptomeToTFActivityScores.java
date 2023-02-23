package org.panda.misc2;

import org.panda.resource.MGI;
import org.panda.utility.FileUtil;
import org.panda.utility.Progress;
import org.panda.utility.statistics.RankedListSignedEnrichment;
import org.panda.utility.statistics.Summary;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class SingleCellTranscriptomeToTFActivityScores
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/home/ozgunbabur/Data/Josh/take3/";
//		convertMouseMatrixToHuman(dir + "normalized_mat.tsv", dir + "normalized_humanized_mat.tsv");
		generate(dir + "normalized_humanized_mat.tsv", dir + "tf_scores.tsv");
//		generate(dir + "simulated_mat.tsv", dir + "simulated_scores.tsv");

		FileUtil.transpose(dir + "tf_scores.tsv", "\t", dir + "tf_scores_t.tsv", "\t", null, null);
	}

	public static void generate(String inFile, String outFile) throws IOException
	{
		String[] cellNames = readCellNames(inFile);
		Map<String, double[]> vMap = readVals(inFile);
		Map<String, double[]> rMap = convertToNormalizedRanks(vMap);
		Map<String, List<String>> rankings = getAllRankings(rMap, cellNames);

		Map<String, Map<String, Boolean>> priors = TFEnrichment.readPriors();
//		Map<String, Map<String, Boolean>> priors = generateSimulationPriors();
		priors = TFEnrichment.convertPriors(priors, priors.values().stream().map(Map::keySet).flatMap(Collection::stream).distinct().collect(Collectors.toList()), 5);

		RankedListSignedEnrichment.writeActivityScoreMatrix(rankings, priors, 3, 1000000, outFile);
	}

	private static Map<String, Map<String, Boolean>> generateSimulationPriors()
	{
		Map<String, Map<String, Boolean>> priors = new HashMap<>();
		Map<String, Boolean> m1 = new HashMap<>();
		for (int i = 0; i < 20; i++)
		{
			m1.put("GT" + i, Boolean.TRUE);
		}
		for (int i = 0; i < 10; i++)
		{
			m1.put("GF" + i, Boolean.TRUE);
		}
		priors.put("TF1", m1);
		Map<String, Boolean> m2 = new HashMap<>();
		for (int i = 0; i < 30; i++)
		{
			m2.put("GF" + (i + 500), Boolean.TRUE);
		}
		priors.put("TF2", m2);
		return priors;
	}

	private static String[] readCellNames(String inFile)
	{
		String line = FileUtil.lines(inFile).findFirst().get();
		line = line.substring(line.indexOf("\t")+1);
		return line.split("\t");
	}

	private static Map<String, double[]> readVals(String inFile)
	{
		Map<String, double[]> valMap = new HashMap<>();

		FileUtil.linesTabbedSkip1(inFile).forEach(t ->
		{
			double[] v = new double[t.length - 1];
			for (int i = 1; i < t.length; i++)
			{
				v[i-1] = Double.parseDouble(t[i]);
			}
			valMap.put(t[0], v);
		});

		return valMap;
	}

	private static Map<String, double[]> convertToZScores(Map<String, double[]> vMap)
	{
		Map<String, double[]> zMap = new HashMap<>();

		vMap.forEach((gene, v) ->
		{
			double mean = Summary.mean(v);

			if (mean == 0) return;

			double sd = Summary.stdev(v);

			double[] z = new double[v.length];
			for (int i = 0; i < v.length; i++)
			{
				z[i] = (v[i] - mean) / sd;
			}
			zMap.put(gene, z);
		});

		return zMap;
	}

	private static Map<String, double[]> convertToNormalizedRanks(Map<String, double[]> vMap)
	{
		Map<String, double[]> rMap = new HashMap<>();

		Progress prg = new Progress(vMap.size(), "Generating normalized ranks");
		vMap.forEach((gene, v) ->
		{
			List<Double> valList = new ArrayList<>();

			for (double value : v)
			{
				if (value > 0) valList.add(value);
			}

			prg.tick();
			if (valList.size() < 10) return;

			valList.sort(Double::compareTo);

			double[] nr = new double[v.length];
			for (int i = 0; i < v.length; i++)
			{
				if (v[i] == 0)
				{
					nr[i] = Double.NaN;
				}
				else
				{
					int fi = valList.indexOf(v[i]);
					int li = valList.lastIndexOf(v[i]);
					double normRank = ((fi + li) / 2D) / valList.size();
					nr[i] = normRank;
				}
			}

			rMap.put(gene, nr);
		});

		return rMap;
	}

	private static List<String> getRankedGenes(Map<String, double[]> zMap, int cellIndex)
	{
		return zMap.keySet().stream().filter(g -> !Double.isNaN(zMap.get(g)[cellIndex]))
			.sorted((g1, g2) -> Double.compare(zMap.get(g2)[cellIndex], zMap.get(g1)[cellIndex]))
			.collect(Collectors.toList());
	}

	private static Map<String, List<String>> getAllRankings(Map<String, double[]> zMap, String[] cellIDs)
	{
		Map<String, List<String>> rankings = new HashMap<>();

		for (int i = 0; i < cellIDs.length; i++)
		{
			rankings.put(cellIDs[i], getRankedGenes(zMap, i));
		}

		return rankings;
	}

	public static void convertMouseMatrixToHuman(String inFile, String outFile) throws IOException
	{
		Map<String, String> lineMap = new HashMap<>();

		Set<String> notUnique = new HashSet<>();

		FileUtil.lines(inFile).skip(1).forEach(l ->
		{
			int tabIndex = l.indexOf("\t");
			String mSym = l.substring(0, tabIndex);
			Set<String> hSyms = MGI.get().getCorrespondingHumanSymbols(mSym);
			if (hSyms.size() == 1)
			{
				String hSym = hSyms.iterator().next();
				if (lineMap.containsKey(hSym))
				{
					notUnique.add(hSym);
				}
				else
				{
					lineMap.put(hSym, l.substring(tabIndex));
				}
			}
		});

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write(FileUtil.readFirstLine(inFile));

		lineMap.keySet().stream().filter(g -> !notUnique.contains(g)).sorted()
			.forEach(g -> FileUtil.lnwrite(g + lineMap.get(g), writer));

		writer.close();
	}

}
