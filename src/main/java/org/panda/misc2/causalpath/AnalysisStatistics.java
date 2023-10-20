package org.panda.misc2.causalpath;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class AnalysisStatistics
{
	public static void main(String[] args)
	{
		String dir = "/home/ozgunbabur/Analyses/Aslan-Thrombin-PAR/";
		String dataFile = dir + "data.csv";
		String[] comp = new String[]{"Resting-vs-PAR1", "Resting-vs-PAR4", "Resting-vs-Thrombin"};

		vennDiagramSignificants(dataFile, "ID", comp, 0.1);

		for (String aCase : comp) countSignificant(dataFile, "ID", aCase);

		dir = dir + "Regular-CP/relax-2aa/";

		for (String aCase : comp) printResultStats(dir + aCase + "/");
	}

	public static void countSignificant(String dataFile, String idCol, String signedPCol)
	{
		System.out.println("dataFile = " + dataFile + "\t" + signedPCol);
		Map<String, String> map = FileUtil.readMap(dataFile, "\t", idCol, signedPCol);
		int rows = map.size();
		long sig = map.keySet().stream().filter(k -> Math.abs(Double.parseDouble(map.get(k))) < 0.1).count();

		System.out.println("rows = " + rows);
		System.out.println("sig = " + sig);
	}

	public static void printResultStats(String dir)
	{
		System.out.println("dir = " + dir);
		String resultFile = dir + "results.txt";
		String[] header = FileUtil.readHeader(resultFile);
		int sInd = ArrayUtil.indexOf(header, "Source data ID");
		int tInd = ArrayUtil.indexOf(header, "Target data ID");

		Set<String> idsContributed = new HashSet<>();
		FileUtil.linesTabbed(resultFile).forEach(t ->
		{
			if (!t[sInd].contains("active")) idsContributed.add(t[sInd]);
			if (!t[tInd].contains("active")) idsContributed.add(t[tInd]);
		});
		System.out.println("Phosphopeptides contributed to at least one edge in the results = " + idsContributed.size());

		Set<String> idsMapped = new HashSet<>();
		FileUtil.linesTabbed(dir + "causative.format").filter(t -> t.length > 3 && t[0].equals("node") && t[2].equals("rppasite")).forEach(t ->
		{
			String id = t[3].substring(0, t[3].indexOf("|"));
			idsMapped.add(id);
		});
		System.out.println("Phosphopeptides mapped on proteins on prior network = " + idsMapped.size());

		Set<String> networkGenes = new HashSet<>();
		Set<String> sourceGenes = new HashSet<>();
		Set<String> targetGenes = new HashSet<>();
		Set<String> loneGenes = new HashSet<>();

		FileUtil.linesTabbed(dir + "causative.sif").forEach(t ->
		{
			if (t.length > 2)
			{
				networkGenes.add(t[0]);
				networkGenes.add(t[2]);
				sourceGenes.add(t[0]);
				targetGenes.add(t[2]);
			}
			else if (t.length == 1)
			{
				loneGenes.add(t[0]);
			}
		});

		loneGenes.removeAll(networkGenes);

		System.out.println("sourceGenes = " + sourceGenes.size());
		System.out.println("targetGenes = " + targetGenes.size());
		System.out.println("networkGenes = " + networkGenes.size());
		System.out.println("loneGenes = " + loneGenes.size());
	}

	public static void vennDiagramSignificants(String dataFile, String idCol, String[] valColNames, double thr)
	{
		String[] header = FileUtil.readHeader(dataFile);
		int idInd = ArrayUtil.indexOf(header, idCol);
		int[] valInds = new int[valColNames.length];
		for (int i = 0; i < valColNames.length; i++)
		{
			valInds[i] = ArrayUtil.indexOf(header, valColNames[i]);
		}

		Set<String>[] groups = new Set[valInds.length];
		for (int i = 0; i < groups.length; i++)
		{
			groups[i] = new HashSet<>();
		}

		Set<String> conflictingIDs = new HashSet<>();

		FileUtil.linesTabbedSkip1(dataFile).forEach(t ->
		{
			String id = t[idInd];

			boolean posit = false;
			boolean negat = false;

			for (int i = 0; i < valInds.length; i++)
			{
				double v = Double.parseDouble(t[valInds[i]]);
				if (Math.abs(v) < thr)
				{
					if (v < 0) negat = true;
					else posit = true;

					groups[i].add(id + (v > 0 ? "+" : "-"));
				}
			}

			if (posit && negat) conflictingIDs.add(id);
		});

		System.out.println("conflictingIDs = " + conflictingIDs);
		System.out.println("-------");

		CollectionUtil.printNameMapping(valColNames);
		CollectionUtil.printVennCounts(groups);

//		System.out.println("-------");
//		for (Set<String> group : groups)
//		{
//			for (String s : group)
//			{
//				System.out.println(s);
//			}
//			System.out.println("\n\n");
//		}
	}

}
