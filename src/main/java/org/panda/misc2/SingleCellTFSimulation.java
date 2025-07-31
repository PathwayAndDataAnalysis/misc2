package org.panda.misc2;

import org.panda.utility.FileUtil;
import org.panda.utility.Progress;
import org.panda.utility.TermCounter;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.Histogram;

import java.io.BufferedWriter;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

public class SingleCellTFSimulation
{
	public static final String DIR = "/home/ozgunbabur/Analyses/CAM/Sim2/";
	static Random rand = new Random();

	public static void main(String[] args)
	{
//		prepareTSVFile();

		reportPerformance("/home/ozgunbabur/Downloads/TF_Activity_Analysis 4.tsv", DIR + "ground_truth.tsv");
	}

	public static void prepareTSVFile()
	{
		int tfCount = 100;
		int geneSize = 20000;
		int cellSize = 10000;
		int tarMin = 5;
		int tarMax = 50;
		int actMin = 100;
		int actMax = 1000;
		double meanMin = 5;
		double meanMax = 20;
		double sdMin = 1;
		double sdMax = 3;
		double shiftAmount = 2;
		double missingValProb = 0.9;

		// Determine TF targets and activity
		Map<String, Map<String, Integer>> tarMap = new HashMap<>();
		Map<String, Map<String, Integer>> celMap = new HashMap<>();

		for (int i = 1; i <= tfCount; i++)
		{
			String tf = "F" + i;

			Map<String, Integer> tars = new HashMap<>();
			int actCnt = randIntUniform(tarMin, tarMax);
			int inhCnt = randIntUniform(tarMin, tarMax);
			Set<Integer> acts = selectRandom(geneSize, actCnt, null);
			Set<Integer> inhs = selectRandom(geneSize, inhCnt, acts);
			acts.forEach(g -> tars.put("G" + g, 1));
			inhs.forEach(g -> tars.put("G" + g, -1));
			tarMap.put(tf, tars);

			Map<String, Integer> cels = new HashMap<>();
			actCnt = randIntUniform(actMin, actMax);
			inhCnt = randIntUniform(actMin, actMax);
			acts = selectRandom(cellSize, actCnt, null);
			inhs = selectRandom(cellSize, inhCnt, acts);
			acts.forEach(c -> cels.put("C" + c, 1));
			inhs.forEach(c -> cels.put("C" + c, -1));
			celMap.put(tf, cels);
		}

		// Save the ground truth
		BufferedWriter writer = FileUtil.newBufferedWriter(DIR + "ground_truth.tsv");
		for (int i = 1; i <= tfCount; i++)
		{
			FileUtil.tab_write("F" + i, writer);
		}
		for (int j = 1; j <= cellSize; j++)
		{
			String cell = "C" + j;
			FileUtil.lnwrite(cell, writer);

			for (int i = 1; i <= tfCount; i++)
			{
				String tf = "F" + i;
				FileUtil.tab_write(celMap.get(tf).getOrDefault(cell, 0), writer);
			}
		}
		FileUtil.closeWriter(writer);

		// Save simulated priors
		BufferedWriter writerP = FileUtil.newBufferedWriter(DIR + "simulated_priors.tsv");
		for (int i = 1; i <= tfCount; i++)
		{
			String tf = "F" + i;

			Map<String, Integer> tars = tarMap.get(tf);
			tars.forEach((tar, dir) -> FileUtil.writeln(tf + "\t" + (dir == 1 ? "up" : "down") + "regulates-expression\t" + tar, writerP));
		}
		FileUtil.closeWriter(writerP);

		// Set gene mean and sd values
		Map<String, Double> meanMap = new HashMap<>();
		Map<String, Double> sdMap = new HashMap<>();
		for (int i = 1; i <= geneSize; i++)
		{
			String gene = "G" + i;
			meanMap.put(gene, randDoubleUniform(meanMin, meanMax));
			sdMap.put(gene, randDoubleUniform(sdMin, sdMax));
		}

		// Generate the RNA data
		DecimalFormat format = new DecimalFormat("#.##");
		writer = FileUtil.newBufferedWriter(DIR + "Simulated_data.tsv");
		for (int i = 1; i <= geneSize; i++)
		{
			FileUtil.tab_write("G" + i, writer);
		}
		Progress prg = new Progress(cellSize, "Generating values");
		for (int j = 1; j <= cellSize; j++)
		{
			prg.tick();
			String cell = "C" + j;
			FileUtil.lnwrite(cell, writer);

			for (int i = 1; i <= geneSize; i++)
			{
				String gene = "G" + i;
				if (rand.nextDouble() < missingValProb)
				{
					FileUtil.tab_write("0", writer);
				} else
				{
					double v = randGaussian(meanMap.get(gene), sdMap.get(gene));
					for (String tf : celMap.keySet())
					{
						int a1 = celMap.get(tf).getOrDefault(cell, 0);
						if (a1 != 0)
						{
							int a2 = tarMap.get(tf).getOrDefault(gene, 0);

							if (a2 != 0)
							{
								int dir = a1 * a2;

								v += dir * (shiftAmount * sdMap.get(gene));
							}
						}
					}
					if (v < 0) v = 0;
					FileUtil.tab_write(format.format(v), writer);
				}
			}
		}
		FileUtil.closeWriter(writer);
	}

	static int randIntUniform(int low, int high)
	{
		return rand.nextInt(high - low + 1) + low;
	}

	static double randDoubleUniform(double low, double high)
	{
		return (rand.nextDouble() * (high - low)) + low;
	}

	static double randGaussian(double mean, double sd)
	{
		return (rand.nextGaussian() * sd) + mean;
	}

	static Set<Integer> selectRandom(int maxNum, int amount, Set<Integer> exclude)
	{
		Set<Integer> set = new HashSet<>();
		while (set.size() < amount)
		{
			int v = randIntUniform(1, maxNum);
			if (exclude == null || !exclude.contains(v)) set.add(v);
		}
		return set;
	}

	static Map<String, Map<String, Double>> readActivitiesDouble(String file)
	{
		Map<String, Map<String, Double>> actMap = new HashMap<>();

		String[] header = FileUtil.readHeader(file);
		for (int i = 1; i < header.length; i++)
		{
			actMap.put(header[i], new HashMap<>());
		}

		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			for (int i = 1; i < t.length; i++)
			{
				double v = Double.NaN;
				if (!t[i].isEmpty()) v = Double.parseDouble(t[i]);
				if (!Double.isNaN(v)) actMap.get(header[i]).put(t[0], v);
			}
		});
		return actMap;
	}

	static Map<String, Map<String, Integer>> readActivitiesInteger(String file)
	{
		Map<String, Map<String, Integer>> actMap = new HashMap<>();

		String[] header = FileUtil.readHeader(file);
		for (int i = 1; i < header.length; i++)
		{
			actMap.put(header[i], new HashMap<>());
		}

		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			for (int i = 1; i < t.length; i++)
			{
//				actMap.get(header[i]).put(t[0], Integer.parseInt(t[i]));
				actMap.get(header[i]).put(t[0], (int) Math.round(Double.parseDouble(t[i])));
			}
		});
		return actMap;
	}

	static Map<String, Map<String, Integer>> discretizeActivities(Map<String, Map<String, Double>> valMap, double fdrThr)
	{
		Map<String, Double> thrMap = new HashMap<>();

		for (String tf : valMap.keySet())
		{
			Map<String, Double> vals = valMap.get(tf);
			Map<String, Double> pvalMap = vals.keySet().stream().collect(Collectors.toMap(c -> c, c -> Math.abs(vals.get(c))));
			double pThr = FDR.getPValueThreshold(pvalMap, null, fdrThr);
			thrMap.put(tf, pThr);
		}

		Map<String, Map<String, Integer>> discMap = new HashMap<>();

		for (String tf : thrMap.keySet())
		{
			Map<String, Integer> map = valMap.get(tf).keySet().stream().collect(Collectors.toMap(c -> c, c ->
				Math.abs(valMap.get(tf).get(c)) <= thrMap.get(tf) ? (valMap.get(tf).get(c) < 0 ? -1 : 1) : 0));
			discMap.put(tf, map);
		}
		return discMap;
	}

	static void reportPerformance(Map<String, Map<String, Integer>> resultMap, Map<String, Map<String, Integer>> truthMap)
	{
		Map<String, Map<String, Integer>> confusionMap = new HashMap<>();
		TermCounter tc = new TermCounter();

		for (String tf : resultMap.keySet())
		{
			Map<String, Integer> cellMap = resultMap.get(tf);
			Map<String, Integer> compMap = truthMap.get(tf);

			Map<String, Integer> counts = new HashMap<>();

			for (String cell : compMap.keySet())
			{
				Integer pred = cellMap.get(cell);
				int truth = compMap.get(cell);

				String type = null;
				if (pred == null)
				{
					type = "NA";
				} else if (pred == 0)
				{
					type = truth == 0 ? "TN" : "FN";
				} else if (pred == 1)
				{
					type = truth == 0 ? "FP" : truth == 1 ? "TP" : "RV";
				} else // (pred == -1)
				{
					type = truth == 0 ? "FP" : truth == -1 ? "TP" : "RV";
				}
				counts.put(type, counts.getOrDefault(type, 0) + 1);
				tc.addTerm(type);
				if (type.equals("RV"))
				{
					System.out.println("cell: " + cell + ", tf: " + tf + ", pred: " + pred + ", truth: " + truth);
				}
			}
			confusionMap.put(tf, counts);
		}
		tc.print();

		System.out.println("Precision = " + (tc.getCountOfTerm("TP") / (double) (tc.getCountOfTerm("TP") + tc.getCountOfTerm("FP"))));
		System.out.println("Recall    = " + (tc.getCountOfTerm("TP") / (double) (tc.getCountOfTerm("TP") + tc.getCountOfTerm("FN"))));


		Map<String, Double> precisionMap = new HashMap<>();
		Map<String, Double> recallMap = new HashMap<>();
		for (String tf : confusionMap.keySet())
		{
			Map<String, Integer> cntMap = confusionMap.get(tf);
			int tp = cntMap.getOrDefault("TP", 0);
			int fp = cntMap.getOrDefault("FP", 0);
			int fn = cntMap.getOrDefault("FN", 0);
			double precision = tp / (double) (tp + fp);
			double recall = tp / (double) (tp + fn);
			precisionMap.put(tf, precision);
			recallMap.put(tf, recall);
		}

		Histogram precisionH = new Histogram(0.1, "Precisions");
		Histogram recallH = new Histogram(0.1, "Recalls");
		precisionH.setBorderAtZero(true);
		recallH.setBorderAtZero(true);
		precisionMap.values().forEach(precisionH::count);
		recallMap.values().forEach(recallH::count);

		System.out.println("\nPrecision histogram:");
		precisionH.print();
		System.out.println("\nRecall histogram:");
		recallH.print();
	}

	static void reportPerformance(String resultFile, String truthFile)
	{
		Map<String, Map<String, Double>> results = readActivitiesDouble(resultFile);
		Map<String, Map<String, Integer>> discRes = discretizeActivities(results, 0.1);
		Map<String, Map<String, Integer>> truth = readActivitiesInteger(truthFile);
		reportPerformance(discRes, truth);
	}
}
