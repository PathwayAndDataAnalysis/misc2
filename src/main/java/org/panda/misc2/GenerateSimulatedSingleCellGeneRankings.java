package org.panda.misc2;

import org.jetbrains.annotations.NotNull;
import org.panda.utility.FileUtil;
import org.panda.utility.Progress;

import java.io.BufferedWriter;
import java.util.*;
import java.util.stream.Collectors;

public class GenerateSimulatedSingleCellGeneRankings
{
	public static final int NUMBER_OF_GENES = 1000;
	public static final int TARGET_MIN = 5;
	public static final int TARGET_MAX = 30;
	public static final double TARGET_ACTIVATED_PROB = 0.5;
	public static final double TF_ACTIVATED_PROB = 0.5;

	public static final int NUMBER_OF_CELLS = 500;
	public static final int CLUSTER_SIZE = 100;
	public static final int TFS_PER_CLUSTER = 10;
	public static final int RANDOM_TF_SIZE = 100;


	public static final String MATRIX_FILE = "/home/ozgunbabur/Documents/Temp/simulated_data_matrix.tsv";
	public static final String PRIORS_FILE = "/home/ozgunbabur/Documents/Temp/simulated_priors.sif";
	public static final String TF_ACT_SIGNS_FILE = "/home/ozgunbabur/Documents/Temp/simulated_tf_act_signs.sif";

	static Random rand = new Random();

	public static void main(String[] args)
	{

		Progress prog = new Progress(NUMBER_OF_CELLS, "Generating ranks");
		Map<String, Map<String, Boolean>> clusterMap = new HashMap<>();
		Map<String, Map<String, Boolean>> tfTargetsMap = new HashMap<>();

		// Generate random TF assignments
		for (int i = 1; i <= RANDOM_TF_SIZE; i++)
		{
			String tfName = "TF-R" + i;
			Map<String, Boolean> targMap = selectTFTargets();
			tfTargetsMap.put(tfName, targMap);
		}

		// Generate 3 clusters, TFs and their targets
		for (int i = 1; i <= 3; i++)
		{
			Map<String, Boolean> clusToTFMap = new HashMap<>();
			String clusterName = "C" + i;
			clusterMap.put(clusterName, clusToTFMap);

			for (int j = 1; j <= TFS_PER_CLUSTER; j++)
			{
				String tfName = "TF-" + clusterName + "-" + j;
				clusToTFMap.put(tfName, Math.random() < TF_ACTIVATED_PROB);

				Map<String, Boolean> targMap = selectTFTargets();
				tfTargetsMap.put(tfName, targMap);
			}
		}

		// Cell ranks
		Map<String, LinkedList<String>> cellRanksMap = new HashMap<>();

		// Generate cluster cells
		for (int i = 1; i <= 3; i++)
		{
			Map<String, Boolean> tfMap = clusterMap.get("C" + i);
			for (int j = 1; j <= CLUSTER_SIZE; j++)
			{
				String name = "C" + i + "-cell-" + j;
				LinkedList<String> list = getRandomizedGenesList();

				for (String tf : tfMap.keySet())
				{
					boolean tfReg = tfMap.get(tf);
					pushGenes(list, tfTargetsMap.get(tf), tfReg, 1);
				}
				cellRanksMap.put(name, list);
				prog.tick();
			}
		}

		// Generate the bridge
		for (int i = 1; i <= CLUSTER_SIZE; i++)
		{
			String name = "B" + "-cell-" + i;
			double leftIntensity = (i - 0.5) / CLUSTER_SIZE;

			LinkedList<String> list = getRandomizedGenesList();

			Map<String, Boolean> c1TFs = clusterMap.get("C1");
			for (String tf : c1TFs.keySet())
			{
				boolean tfReg = c1TFs.get(tf);
				pushGenes(list, tfTargetsMap.get(tf), tfReg, leftIntensity);
			}

			Map<String, Boolean> c2TFs = clusterMap.get("C2");
			for (String tf : c2TFs.keySet())
			{
				boolean tfReg = c2TFs.get(tf);
				pushGenes(list, tfTargetsMap.get(tf), tfReg, 1 - leftIntensity);
			}

			cellRanksMap.put(name, list);
			prog.tick();
		}

		// Generate random cells

		for (int i = 1; i <= NUMBER_OF_CELLS - (4 * CLUSTER_SIZE); i++)
		{
			String name = "R" + "-cell-" + i;
			LinkedList<String> list = getRandomizedGenesList();
			cellRanksMap.put(name, list);
			prog.tick();
		}

		// Write the data matrix
		List<String> cells = cellRanksMap.keySet().stream().sorted().collect(Collectors.toList());
		Map<String, Map<String, Double>> geneToCellToRank = new HashMap<>();
		for (String cell : cellRanksMap.keySet())
		{
			LinkedList<String> genes = cellRanksMap.get(cell);
			int ind = 0;
			for (Iterator<String> it = genes.iterator(); it.hasNext(); )
			{
				String gene = it.next();
				double rank = 1 - (((ind++) + 0.5) / genes.size());
				if (!geneToCellToRank.containsKey(gene)) geneToCellToRank.put(gene, new HashMap<>());
				geneToCellToRank.get(gene).put(cell, rank);
			}
		}

		prog = new Progress(NUMBER_OF_GENES, "Writing data matrix");

		BufferedWriter writer = FileUtil.newBufferedWriter(MATRIX_FILE);
		for (String cell : cells)
		{
			FileUtil.tab_write(cell, writer);
		}
		for (int i = 1; i <= NUMBER_OF_GENES; i++)
		{
			prog.tick();
			String gene = "G" + i;
			FileUtil.lnwrite(gene, writer);

			Map<String, Double> cellToRank = geneToCellToRank.get(gene);

			for (String cell : cells)
			{
				double rank = cellToRank.get(cell);
				FileUtil.tab_write(rank, writer);
			}
		}
		FileUtil.closeWriter(writer);

		// Write the priors file
		writer = FileUtil.newBufferedWriter(PRIORS_FILE);
		for (String tf : tfTargetsMap.keySet())
		{
			Map<String, Boolean> targMap = tfTargetsMap.get(tf);
			for (String targ : targMap.keySet())
			{
				boolean sign = targMap.get(targ);
				FileUtil.writeln(tf + "\t" + (sign ? "up" : "down") + "regulates-expression\t" + targ, writer);
			}
		}
		FileUtil.closeWriter(writer);

		writer = FileUtil.newBufferedWriter(TF_ACT_SIGNS_FILE);
		FileUtil.write("Cluster\tTF name\tDirection", writer);

		// Write TF directions for each cluster
		for (String clusterName : clusterMap.keySet())
		{
			Map<String, Boolean> tfActMap = clusterMap.get(clusterName);
			for (String tf : tfActMap.keySet())
			{
				FileUtil.lnwrite(clusterName + "\t" + tf + "\t" + (tfActMap.get(tf) ? "+" : "-"), writer);
			}
		}
		FileUtil.closeWriter(writer);
	}

	@NotNull
	private static Map<String, Boolean> selectTFTargets()
	{
		int targetSize = rand.nextInt(TARGET_MAX - TARGET_MIN + 1) + TARGET_MIN;

		Set<Integer> targInds = new HashSet<>();
		while (targInds.size() < targetSize)
		{
			int ind = rand.nextInt(NUMBER_OF_GENES) + 1;
			targInds.add(ind);
		}

		Map<String, Boolean> targMap = new HashMap<>();
		for (Integer ind : targInds)
		{
			targMap.put("G" + ind, rand.nextDouble() < TARGET_ACTIVATED_PROB);
		}
		return targMap;
	}

	@NotNull
	private static LinkedList<String> getRandomizedGenesList()
	{
		LinkedList<String> list = new LinkedList<>();
		for (int k = 1; k <= NUMBER_OF_GENES; k++) list.add("G" + k);
		Collections.shuffle(list);
		return list;
	}

	static void pushGenes(LinkedList<String> genes, Map<String, Boolean> targetsMap, boolean tfReg, double intensity)
	{
		for (String gene : targetsMap.keySet())
		{
			boolean reg = targetsMap.get(gene);
			int ind = genes.indexOf(gene);

			if (ind < 0)
			{
				System.out.println("What's going on??? ind = " + ind);
			}

			int newInd = -1;

			if (tfReg ^ reg) // push to the end
			{
				int range = (int) ((NUMBER_OF_GENES - ind) * intensity);
				if (range > 0)
				{
					int move = rand.nextInt(range);
					newInd = ind + move - 1;
				}
			}
			else // push to front
			{
				int bound = (int) (ind * intensity);
				if (bound > 0)
				{
					int move = rand.nextInt(bound);
					newInd = ind - move;
				}
			}

			if (newInd >= 0)
			{
				genes.remove(ind);
				genes.add(newInd, gene);
			}
		}
	}
}
