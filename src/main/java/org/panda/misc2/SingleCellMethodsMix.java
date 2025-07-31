package org.panda.misc2;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.panda.utility.*;
import org.panda.utility.statistics.Binomial;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.StouffersCombinedProbability;
import org.panda.utility.statistics.TTest;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class SingleCellMethodsMix
{
	public static final String BASE = "/home/ozgunbabur/Data/scPerturb/";

//	public static final String STUDY = "ReplogleWeissman2022_K562_essential";
//	public static final String STUDY = "ReplogleWeissman2022_rpe1";
	public static final String STUDY = "FrangiehIzar2021_RNA";

//	public static final String[] TRANSFER_STUDIES = new String[]{"ReplogleWeissman2022_rpe1", "FrangiehIzar2021_RNA"};
//	public static final String[] TRANSFER_STUDIES = new String[]{"ReplogleWeissman2022_K562_essential", "FrangiehIzar2021_RNA"};
	public static final String[] TRANSFER_STUDIES = new String[]{"ReplogleWeissman2022_K562_essential", "ReplogleWeissman2022_rpe1"};

	public static final String ZSCORES = BASE + "z-scores-v2/" + STUDY + "_zscores.csv";
	public static final String BIN_FILE = BASE + "binary/" + STUDY + ".bin";
	public static final String META = BASE + "cell-metadata/" + STUDY + "-metadata.csv";

//	public static final String TF = "MAX";
	public static final String TF = "MYC";
//	public static final String TF = "CCND1";
//	public static final String TF = "DNMT1";

	public static int MIN_N = 100;
	public static int MIN_K = 5;
	public static final double FDR_THR = 0.1;
	public static int SIGNATURE_SIZE = 100;


	public static void main(String[] args)
	{
		System.out.println("STUDY = " + STUDY);
//		System.out.println("TF = " + TF);

//		Map<String, double[]> dataset = loadZScores(ZSCORES);
//		saveDataset(dataset, BIN_FILE);

//		writeTTests();
//		classify();
		testAllTFs();

	}

	private static void testAllTFs()
	{
		Map<String, double[]> dataset = readBinaryDataset(BIN_FILE);//loadZScores(ZSCORES);

		int cnum = dataset.values().iterator().next().length;
		System.out.println("cell = " + cnum);
		System.out.println("genes = " + dataset.size());
		Map<String, Map<String, Integer>> priors = loadPriors();

		String line = FileUtil.readFirstLine(ZSCORES);
		line = line.substring(line.indexOf("\t")+1);
		String[] cells = line.split("\t");

		for (String tf : priors.keySet())
		{
			Map<String, Integer> targMap = priors.get(tf);
			if (targMap.size() >= 5)
			{
				Set<String> cellNames = getSubsetColNames(META, tf);

				if (!cellNames.isEmpty())
				{
					System.out.println("tf = " + tf);
					int[] cInds = getInds(cells, cellNames);

					Map<String, Double> pMap = new HashMap<>();

					doPlainRanksumForInactivation(dataset, cells, cInds, targMap, pMap);
					double p = combinePValues(pMap);
					System.out.println("p = " + p);
				}
			}
		}
	}

	private static void writeTTests()
	{
		Map<String, double[]> dataset = readBinaryDataset(BIN_FILE);//loadZScores(ZSCORES);
		System.out.println("dataset.size() = " + dataset.size());
		Set<String> controls = getSubsetColNames(META, "control");
		System.out.println("controls.size() = " + controls.size());
		Set<String> tests = getSubsetColNames(META, TF);
		System.out.println("tests.size() = " + tests.size());

		String line = FileUtil.readFirstLine(ZSCORES);
		line = line.substring(line.indexOf("\t")+1);
		String[] cells = line.split("\t");

		Map<String, Tuple> tTests = getTTests(dataset, cells, controls, tests);
		System.out.println("tTests.size() = " + tTests.size());

		BufferedWriter writer = FileUtil.newBufferedWriter(getTTestFile(STUDY, TF));
		FileUtil.write("Gene\tt-val\tp-val", writer);
		tTests.forEach((s, tuple) -> FileUtil.lnwrite(s + "\t" + -tuple.v + "\t" + tuple.p, writer));
		FileUtil.closeWriter(writer);
	}

	private static void classify()
	{
		Map<String, double[]> dataset = readBinaryDataset(BIN_FILE);//loadZScores(ZSCORES);

		int cnum = dataset.values().iterator().next().length;
		System.out.println("cnum = " + cnum);

		System.out.println("dataset.size() = " + dataset.size());
		Set<String> controls = getSubsetColNames(META, "control");
		System.out.println("controls.size() = " + controls.size());
		Set<String> tests = getSubsetColNames(META, TF);
		System.out.println("tests.size() = " + tests.size());
		Set<String> subsetCells = new HashSet<>(controls);
		subsetCells.addAll(tests);

		String line = FileUtil.readFirstLine(ZSCORES);
		line = line.substring(line.indexOf("\t")+1);
		String[] cells = line.split("\t");
		int[] cInds = getInds(cells, subsetCells);

		Map<String, Integer> signMap = new HashMap<>();
		Map<String, Double> pMap = new HashMap<>();

		System.out.println("\nPlain ranksum");
		Map<String, Map<String, Integer>> priorMap = loadPriors();
		Map<String, Integer> tars = priorMap.get(TF);
		doPlainRanksum(dataset, cells, cInds, tars, signMap, pMap);
		printClassificationResults(controls, tests, signMap, pMap);

		System.out.println("\nWeighted ranksum using self-weights");
		Map<String, Tuple> tTests = loadTTests(getTTestFile(STUDY, TF));
		Map<String, Double> preWeights = getPreWeights(tTests, tars);
		doWeightedRanksum(dataset, cells, cInds, preWeights, signMap, pMap);
		printClassificationResults(controls, tests, signMap, pMap);

		System.out.println("\nWeighted ranksum using transfer weights");
		for (String transferStudy : TRANSFER_STUDIES)
		{
			System.out.println("\ntransferStudy = " + transferStudy);
			Map<String, Tuple> transTTests = loadTTests(getTTestFile(transferStudy, TF));
			preWeights = getPreWeights(transTTests, tars);
			doWeightedRanksum(dataset, cells, cInds, preWeights, signMap, pMap);
			printClassificationResults(controls, tests, signMap, pMap);
		}

		System.out.println("\nSignature ranksum using self signature");
		preWeights = getPreWeights(tTests, SIGNATURE_SIZE);
		doWeightedRanksum(dataset, cells, cInds, preWeights, signMap, pMap);
		printClassificationResults(controls, tests, signMap, pMap);

		System.out.println("\nSignature ranksum using transfer signature");
		for (String transferStudy : TRANSFER_STUDIES)
		{
			System.out.println("\ntransferStudy = " + transferStudy);
			Map<String, Tuple> transTTests = loadTTests(getTTestFile(transferStudy, TF));
			preWeights = getPreWeights(transTTests, SIGNATURE_SIZE);
			doWeightedRanksum(dataset, cells, cInds, preWeights, signMap, pMap);
			printClassificationResults(controls, tests, signMap, pMap);
		}
	}

	private static void printClassificationResults(Set<String> cells, Map<String, Integer> signMap, Map<String, Double> pMap)
	{
		List<String> sigCells = FDR.select(pMap, null, FDR_THR);
		int[] cnt = new int[]{0, 0, 0};

		for (String cell : cells)
		{
			if (sigCells.contains(cell))
			{
				if (signMap.get(cell) > 0) cnt[0]++;
				else cnt[2]++;
			}
			else cnt[1]++;
		}
		System.out.println(cnt[0] + "\t" + cnt[1] + "\t" + cnt[2]);
	}

	private static void analyzeSignDistribution(Map<String, Integer> signMap)
	{
		int pos = (int) signMap.keySet().stream().filter(n -> signMap.get(n) > 0).count();
		int neg = (int) signMap.keySet().stream().filter(n -> signMap.get(n) < 0).count();
		double pval = Binomial.getPval(pos, neg);
		System.out.println("<" + pos + " | " + neg + ">");
		System.out.println("pval = " + pval);
	}

	private static double combinePValues(Map<String, Double> pvals)
	{
		double[] p = new double[pvals.size()];
		int i = 0;
		for (String cell : pvals.keySet())
		{
			p[i++] = pvals.get(cell);
		}
		return StouffersCombinedProbability.combineP(p);
	}

	private static void printClassificationResults(Set<String> controls, Set<String> tests, Map<String, Integer> signMap, Map<String, Double> pMap)
	{
		List<String> sigCells = FDR.select(pMap, null, FDR_THR);

		int[][] cnt = new int[][]{{0, 0, 0},{0, 0, 0}};

		for (String ctrl : controls)
		{
			if (sigCells.contains(ctrl))
			{
				if (signMap.get(ctrl) > 0) cnt[0][0]++;
				else cnt[0][2]++;
			}
			else cnt[0][1]++;
		}
		for (String test : tests)
		{
			if (sigCells.contains(test))
			{
				if (signMap.get(test) > 0) cnt[1][0]++;
				else cnt[1][2]++;
			}
			else cnt[1][1]++;
		}

		System.out.println(cnt[0][0] + "\t" + cnt[0][1] + "\t" + cnt[0][2]);
		System.out.println(cnt[1][0] + "\t" + cnt[1][1] + "\t" + cnt[1][2]);

		double[][] rats = new double[2][3];
		for (int i = 0; i < 3; i++)
		{
			rats[0][i] = cnt[0][i] / (double) (ArrayUtil.sum(cnt[0]));
		}
		for (int i = 0; i < 3; i++)
		{
			rats[1][i] = cnt[1][i] / (double) (ArrayUtil.sum(cnt[1]));
		}
		System.out.println(rats[0][0] + "\t" + rats[0][1] + "\t" + rats[0][2]);
		System.out.println(rats[1][0] + "\t" + rats[1][1] + "\t" + rats[1][2]);
	}

	private static void doPlainRanksum(Map<String, double[]> dataset, String[] cells, int[] cInds, Map<String, Integer> tars, Map<String, Integer> signMap, Map<String, Double> pMap)
	{
		for (int ind : cInds)
		{
			List<String> rankedGenes = getRankedGenes(dataset, ind);
			int n = rankedGenes.size();
			Map<String, Double> normalizedRanks = getNormalizedRanks(rankedGenes);
			int k = (int) tars.keySet().stream().filter(normalizedRanks::containsKey).count();
//			if (n < MIN_N || k < MIN_K) continue;
			if (k == 0) continue;
			double meanRank = getMeanRank(normalizedRanks, tars);
			signMap.put(cells[ind], meanRank <= 0.5 ? 1 : -1);
			double pval = getSignificance(meanRank, n, k);
			pMap.put(cells[ind], pval);
		}
	}

	private static void doPlainRanksumForInactivation(Map<String, double[]> dataset, String[] cells, int[] cInds, Map<String, Integer> tars, Map<String, Double> pMap)
	{
		for (int ind : cInds)
		{
			List<String> rankedGenes = getRankedGenes(dataset, ind);
			int n = rankedGenes.size();
			Map<String, Double> normalizedRanks = getNormalizedRanks(rankedGenes);
			int k = (int) tars.keySet().stream().filter(normalizedRanks::containsKey).count();
			if (k == 0) continue;
			double meanRank = 1 - getMeanRank(normalizedRanks, tars);
			double pval = getSignificanceOneTailed(meanRank, n, k);
			pMap.put(cells[ind], pval);
		}
	}

	private static void doWeightedRanksum(Map<String, double[]> dataset, String[] cells, int[] cInds, Map<String, Double> tarPreWeights, Map<String, Integer> signMap, Map<String, Double> pMap)
	{
		for (int ind : cInds)
		{
			List<String> rankedGenes = getRankedGenes(dataset, ind);
			int n = rankedGenes.size();
			Map<String, Double> normalizedRanks = getNormalizedRanks(rankedGenes);
			Map<String, Double> weights = adjustWeightsToRanking(normalizedRanks, tarPreWeights);

			if (n < MIN_N || weights.size() < MIN_K) continue;

			double meanRank = getWeightedMeanRank(normalizedRanks, weights);
			signMap.put(cells[ind], meanRank <= 0.5 ? 1 : -1);
			double pval = getSignificance(meanRank, n, weights);
			pMap.put(cells[ind], pval);
		}
	}

	public static Map<String, double[]> loadZScores(String file)
	{
		Map<String, double[]> dataset = new HashMap<>();
		String[] header = FileUtil.readHeader(file);
		int cellsSize = header.length - 1;

		int lineNum = (int) FileUtil.findNumberOfLines(file);
		System.out.println("lineNum = " + lineNum);

		Progress prg = new Progress(lineNum - 1, "Reading file");
		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			prg.tick();
			double[] row = new double[cellsSize];
			dataset.put(t[0], row);

			for (int i = 1; i < t.length; i++)
			{
				row[i-1] = t[i].isEmpty() ? Double.NaN : Double.parseDouble(t[i]);
			}
		});
		return dataset;
	}


	public static List<String> getRankedGenes(Map<String, double[]> dataset, int cell)
	{
		List<String> genes = new ArrayList<>(dataset.keySet());
		Set<String> dropped = new HashSet<>();
		genes.forEach(gene ->
		{
			if (Double.isNaN(dataset.get(gene)[cell])) dropped.add(gene);
		});
		genes.removeAll(dropped);
		genes.sort(Comparator.comparing(g -> -dataset.get(g)[cell]));
		return genes;
	}

	public static Map<String, Double> getNormalizedRanks(List<String> genes)
	{
		Map<String, Double> normRanks = new HashMap<>();
		int n = genes.size();
		for (int i = 0; i < genes.size(); i++)
		{
			normRanks.put(genes.get(i), (i + 0.5) / n);
		}
		return normRanks;
	}

	public static double getMeanRank(Map<String, Double> normRanks, Map<String, Integer> tars)
	{
		int k = 0;
		double sum = 0;

		for (String targ : tars.keySet())
		{
			if (normRanks.containsKey(targ))
			{
				sum += tars.get(targ) > 0 ? normRanks.get(targ) : 1 - normRanks.get(targ);
				k++;
			}
		}
		return sum / k;
	}

	public static Map<String, Double> adjustWeightsToRanking(Map<String, Double> ranking, Map<String, Double> preWeigts)
	{
		Set<String> genes = ranking.keySet();
		double sum = 0;
		Set<String> existing = new HashSet<>();
		for (String gene : preWeigts.keySet())
		{
			if (genes.contains(gene))
			{
				sum += Math.abs(preWeigts.get(gene));
				existing.add(gene);
			}
		}
		Map<String, Double> adjustedW = new HashMap<>();
		for (String gene : existing)
		{
			adjustedW.put(gene, preWeigts.get(gene) / sum);
		}
		return adjustedW;
	}

	public static double getWeightedMeanRank(Map<String, Double> normRanks, Map<String, Double> weights)
	{
		double mr = 0;

		for (String targ : weights.keySet())
		{
			double w = weights.get(targ);
			double rank = normRanks.get(targ);
			mr += (w > 0 ? rank : 1 - rank) * Math.abs(w);
		}
		return mr;
	}

	public static double getSignificance(double meanRank, double n, double k)
	{
		if (meanRank > 0.5) meanRank = 1 - meanRank;

		double sd = Math.sqrt(((n + 1) * (n-k)) / (12 * n * n * k));

		NormalDistribution nd = new NormalDistribution(0.5, sd);
		return 2 * nd.cumulativeProbability(meanRank);
	}

	public static double getSignificanceOneTailed(double meanRank, double n, double k)
	{
		double sd = Math.sqrt(((n + 1) * (n-k)) / (12 * n * n * k));

		NormalDistribution nd = new NormalDistribution(0.5, sd);
		return nd.cumulativeProbability(meanRank);
	}

	public static double getSignificance(double meanRank, double n, Map<String, Double> weights)
	{
		if (meanRank > 0.5) meanRank = 1 - meanRank;
		double s = 0;
		for (String gene : weights.keySet())
		{
			double w = weights.get(gene);
			s += w * w;
		}

		double sd = Math.sqrt(((n + 1) * ((n * s) - 1)) / (12 * n * n));

		NormalDistribution nd = new NormalDistribution(0.5, sd);
		return 2 * nd.cumulativeProbability(meanRank);
	}

	public static Tuple doTTest(Map<String, Double> rowMap, Set<String> controls, Set<String> tests)
	{
		List<Double> cList = new ArrayList<>();
		List<Double> tList = new ArrayList<>();
		for (String control : controls)
		{
			if (rowMap.containsKey(control)) cList.add(rowMap.get(control));
		}
		for (String test : tests)
		{
			if (rowMap.containsKey(test)) tList.add(rowMap.get(test));
		}
		double[] ctrlVals = new double[cList.size()];
		for (int i = 0; i < cList.size(); i++)
		{
			ctrlVals[i] = cList.get(i);
		}
		double[] testVals = new double[tList.size()];
		for (int i = 0; i < tList.size(); i++)
		{
			testVals[i] = tList.get(i);
		}
		return new Tuple(TestUtils.t(ctrlVals, testVals), TestUtils.tTest(ctrlVals, testVals));
	}

	public static Map<String, Tuple> loadTTests(String filename)
	{
		Map<String, Tuple> map = new HashMap<>();
		FileUtil.linesTabbedSkip1(filename).forEach(t ->
			map.put(t[0], new Tuple(Double.parseDouble(t[1]), Double.parseDouble(t[2]))));
		return map;
	}

	public static Map<String, Double> getPreWeights(Map<String, Tuple> ttests, int size)
	{
//		size = Math.min(size, findSignificant(ttests, FDR_THR).size());
		List<String> list = new ArrayList<>(ttests.keySet());
		list.sort(Comparator.comparing(g -> ttests.get(g).p));
		Map<String, Double> preW = new HashMap<>();
		for (int i = 0; i < size; i++)
		{
			String gene = list.get(i);
			Tuple tup = ttests.get(gene);
			preW.put(gene, Math.signum(tup.v) * -Math.log(tup.p));
		}
		System.out.println("preW.size() = " + preW.size());
		return preW;
	}

	public static Map<String, Double> getPreWeights(Map<String, Tuple> ttests, Map<String, Integer> tars)
	{
		Map<String, Double> preW = new HashMap<>();
		for (String tar : tars.keySet())
		{
			if (ttests.containsKey(tar))
			{
				Tuple tup = ttests.get(tar);
				if (tup.v * tars.get(tar) > 0)
				{
					preW.put(tar, Math.signum(tup.v) * -Math.log(tup.p));
				}
			}
		}
		return preW;
	}

	public static Map<String, Tuple> getSignificantGenes(Map<String, Map<String, Double>> dataset, Set<String> controls, Set<String> tests, double fdrThr)
	{
		Map<String, Tuple> map = new HashMap<>();
		Progress prg = new Progress(dataset.size(), "Conducting t-tests");
		for (String gene : dataset.keySet())
		{
			Tuple tuple = doTTest(dataset.get(gene), controls, tests);
			map.put(gene, tuple);
			prg.tick();
		}
		Map<Tuple, Double> pvals = map.values().stream().collect(Collectors.toMap(t -> t, t -> t.p));
		double pThr = FDR.getPValueThreshold(pvals, null, fdrThr);
		List<String> toRemove = new ArrayList<>();
		map.forEach((gene, tuple) -> {
			if (tuple.p > pThr) toRemove.add(gene);
		});
		toRemove.forEach(map::remove);
		return map;
	}

	public static Set<String> getSubsetColNames(String metadataFile, String perturb)
	{
		String[] header = FileUtil.readHeader(metadataFile);
		int pertInd = ArrayUtil.indexOf(header, "perturbation");

		return FileUtil.linesTabbedSkip1(metadataFile)
			.filter(t -> t[pertInd].equals(perturb))
			.map(t -> t[0]).collect(Collectors.toSet());
	}

	public static Map<String, Tuple> getTTests(Map<String, double[]> dataset, String[] cells, Set<String> controls, Set<String> tests)
	{
		Map<String, Tuple> tMap = new HashMap<>();

		int[] cInds = getInds(cells, controls);
		int[] tInds = getInds(cells, tests);

		Progress prg = new Progress(dataset.size(), "Testing");
		dataset.keySet().forEach(gene ->
		{
			prg.tick();
			double[] row = dataset.get(gene);

			List<Double> ctrlList = new ArrayList<>();
			List<Double> testList = new ArrayList<>();

			for (int i = 0; i < cInds.length; i++)
			{
				if (!Double.isNaN(row[cInds[i]])) ctrlList.add(row[cInds[i]]);
			}
			for (int i = 0; i < tInds.length; i++)
			{
				if (!Double.isNaN(row[tInds[i]])) testList.add(row[tInds[i]]);
			}

			if (ctrlList.size() < 3 || testList.size() < 3) return;

			double[] ctrlVals = new double[ctrlList.size()];
			double[] testVals = new double[testList.size()];
			for (int i = 0; i < ctrlVals.length; i++)
			{
				ctrlVals[i] = ctrlList.get(i);
			}
			for (int i = 0; i < testVals.length; i++)
			{
				testVals[i] = testList.get(i);
			}

			tMap.put(gene, TTest.test(ctrlVals, testVals));
		});
		return tMap;
	}

	public static int[] getInds(String[] array, Set<String> select)
	{
		int[] inds = new int[select.size()];
		int i = 0;
		for (int j = 0; j < array.length; j++)
		{
			if (select.contains(array[j])) inds[i++] = j;
		}
		return inds;
	}

	public static Map<String, Tuple> findSignificant(Map<String, Tuple> tMap, double fdrThr)
	{
		Map<String, Double> pMap = new HashMap<>();
		tMap.keySet().forEach(gene -> pMap.put(gene, tMap.get(gene).p));
		List<String> select = FDR.select(pMap, null, fdrThr);
		Map<String, Tuple> sMap = new HashMap<>();
		select.forEach(gene -> sMap.put(gene, tMap.get(gene)));
		return sMap;
	}

	public static Map<String, Double> generateWeights(Map<String, Tuple> tTestResults, Set<String> select)
	{
		Map<String, Double> rawW = new HashMap<>();
		double sum = 0;
		for (String gene : select)
		{
			if (tTestResults.containsKey(gene))
			{
				double rw = -Math.log(tTestResults.get(gene).p);
				sum += rw;
				rawW.put(gene, rw);
			}
		}
		Map<String, Double> weights = new HashMap<>();
		for (String gene : rawW.keySet())
		{
			weights.put(gene, rawW.get(gene) / sum);
		}
		return weights;
	}

	public static Map<String, Map<String, Integer>> loadPriors()
	{
		Map<String, Map<String, Integer>> priors = new HashMap<>();
		FileUtil.linesTabbed("/home/ozgunbabur/Data/causal-priors.txt")
			.filter(t -> t[1].equals("upregulates-expression") || t[1].equals("downregulates-expression")).forEach(t ->
			{
				if (!priors.containsKey(t[0])) priors.put(t[0], new HashMap<>());
				priors.get(t[0]).put(t[2], t[1].startsWith("u") ? 1 : -1);
			});
		return priors;
	}

	public static void saveDataset(Map<String, double[]> dataset, String filename)
	{
		try
		{
			ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(filename));
			oos.writeObject(dataset);
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	public static Map<String, double[]> readBinaryDataset(String filename)
	{
		try
		{
			Kronometre kro = new Kronometre();
			ObjectInputStream ois = new ObjectInputStream(new FileInputStream(filename));
			Map<String, double[]> dataset = (Map<String, double[]>) ois.readObject();
			kro.stop();
			kro.print();
			return dataset;
		}
		catch (IOException | ClassNotFoundException e)
		{
			e.printStackTrace();
			return null;
		}
	}

	public static String getTTestFile(String study, String tf)
	{
		return BASE + "ttests/" + study + "-" + tf + "-ttests.tsv";
	}
}
