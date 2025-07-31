package org.panda.misc2.siffile;

import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class DownstreamPathExtractor
{
	public static void main(String[] args) throws IOException
	{
//		String dir = "/home/ozgunbabur/Analyses/Platelet-Blood-paper/merged/";
//		extract(dir + "merged-relax1aa.sif", dir + "GP6-downstream.sif", "GP6");

		String dir = "/home/ozgunbabur/Analyses/CPTAC-LSCC/v3/altered-subset-vs-normals/KEAP1-NFE2L2-CUL3-altered/";

//		extractNaive(dir + "causative.sif", dir + "NFE2L2-downstream.sif", 2, "NFE2L2");

		Map<String, Integer> src = new HashMap<>();
		src.put("NFE2L2", 1);
		extractCausalFromResultsFile(dir, dir + "NFE2L2-downstream-2.sif", 2, src);
	}

	public static void extractNaive(String inFile, String outFile, int limit, String... sources) throws IOException
	{
		Set<String> lines = FileUtil.lines(inFile).filter(l -> l.split("\t").length > 2).collect(Collectors.toSet());

		Set<String> subsetLines = new HashSet<>();

		Set<String> protsToFollow = new HashSet<>();
		Set<String> nextProtsToFollow = new HashSet<>(Arrays.asList(sources));
		Set<String> visited = new HashSet<>();

		int iter = 0;

		while (!nextProtsToFollow.isEmpty() && iter < limit)
		{
			protsToFollow.addAll(nextProtsToFollow);
			nextProtsToFollow.clear();

			System.out.println("protsToFollow = " + protsToFollow);

			for (String line : lines)
			{
				String[] t = line.split("\t");
				if (protsToFollow.contains(t[0]))
				{
					if (!visited.contains(t[2]) && !protsToFollow.contains(t[2]))
					{
						subsetLines.add(line);
						nextProtsToFollow.add(t[2]);
					}
				}
			}

			visited.addAll(protsToFollow);
			protsToFollow.clear();

			iter++;
		}

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		subsetLines.forEach(l -> FileUtil.writeln(l, writer));
		FileUtil.closeWriter(writer);

		if (inFile.endsWith(".sif") && outFile.endsWith(".sif"))
		{
			String formatF = inFile.substring(0, inFile.length() - 3) + "format";
			if (FileUtil.exists(formatF))
			{
				String outFormat = outFile.substring(0, outFile.length() - 3) + "format";
				FileUtil.copyFile(formatF, outFormat);
			}
		}
	}

	private static void extractCausalFromResultsFile(String dir, String outSif, int limit, Map<String, Integer> sourceWithAct) throws IOException
	{
		String resultsFile = dir + "results.txt";
		String[] header = FileUtil.readHeader(resultsFile);
		int sourceInd = ArrayUtil.indexOf(header, "Source");
		int relInd = ArrayUtil.indexOf(header, "Relation");
		int targetInd = ArrayUtil.indexOf(header, "Target");
		int srcEffectInd = ArrayUtil.indexOf(header, "Source site effect");
		int tgtEffectInd = ArrayUtil.indexOf(header, "Target site effect");
		int srcChangeInd = ArrayUtil.indexOf(header, "Source change");
		int tgtChangeInd = ArrayUtil.indexOf(header, "Target change");
		int srcTypeInd = ArrayUtil.indexOf(header, "Source data type");
		int tgtTypeInd = ArrayUtil.indexOf(header, "Target data type");
		int srcID = ArrayUtil.indexOf(header, "Source data ID");
		int tgtID = ArrayUtil.indexOf(header, "Target data ID");

		Map<String, String> sifMap = FileUtil.lines(dir + "causative.sif")
			.filter(l -> l.split("\t").length > 2)
			.collect(Collectors.toMap(l -> {
				String[] t = l.split("\t");
				return t[0] + "\t" + t[1] + "\t" + t[2];
				}, l -> l));

		Set<String[]> resRows = FileUtil.linesTabbed(resultsFile).skip(1).collect(Collectors.toSet());

		Set<String> subsetLines = new HashSet<>();
		Set<String> ids = new HashSet<>();
		Set<String> syms = new HashSet<>();

		Map<String, Integer> protsToFollow = new HashMap<>();
		Map<String, Integer> nextProtsToFollow = new HashMap<>();
		nextProtsToFollow.putAll(sourceWithAct);
		Set<String> visited = new HashSet<>();

		int iter = 0;

		while (!nextProtsToFollow.isEmpty() && iter < limit)
		{
			protsToFollow.putAll(nextProtsToFollow);
			nextProtsToFollow.clear();

			System.out.println("protsToFollow = " + protsToFollow);

			for (String[] t : resRows)
			{
				if (protsToFollow.containsKey(t[0]))
				{
					if (getSourceActivity(t, srcTypeInd, srcEffectInd, srcChangeInd) == protsToFollow.get(t[0]) &&
						!visited.contains(t[2])) //&& !protsToFollow.containsKey(t[2]))
					{
						String key = t[sourceInd] + "\t" + t[relInd] + "\t" + t[targetInd];
						subsetLines.add(sifMap.get(key));

						ids.add(t[srcID]);
						ids.add(t[tgtID]);
						syms.add(t[sourceInd]);
						syms.add(t[targetInd]);

						int targetAct = getTargetActivity(t, tgtTypeInd, tgtEffectInd, tgtChangeInd);
						if (targetAct != 0)
						{
							nextProtsToFollow.put(t[2], targetAct);
						}
					}
				}
			}
			visited.addAll(protsToFollow.keySet());
			protsToFollow.clear();

			iter++;
		}

		BufferedWriter writer = FileUtil.newBufferedWriter(outSif);
		subsetLines.forEach(l -> FileUtil.writeln(l, writer));
		FileUtil.closeWriter(writer);

		String formatF = dir + "causative.format";
		String outFormat = outSif.substring(0, outSif.length() - 3) + "format";
		BufferedWriter fWriter = FileUtil.newBufferedWriter(outFormat);

		// Copy the subset format, exclude unknown-effect modif-proteomics if they are not in causal reasoning.
		FileUtil.lines(formatF).forEach(l ->
		{
			String[] t = l.split("\t");

			// condition to skip
			if (t[0].equals("node") && (
				(!t[1].equals("all-nodes") && !syms.contains(t[1])) ||
				(t[2].equals("rppasite") && !ids.contains(t[3].substring(0, t[3].indexOf("|"))) &&
					(t[3].contains("|p|") || t[3].contains("|a|") || t[3].contains("|m|") || t[3].contains("|u|")) && t[3].contains("|50 50 50|"))
			)) return;

			FileUtil.writeln(l, fWriter);
		});

		FileUtil.closeWriter(fWriter);
	}

	private static int getSourceActivity(String[] t, int srcTypeInd, int srcEffectInd, int srcChangeInd)
	{
		int effect = t[srcTypeInd].endsWith("protein") && t[srcTypeInd].length() > 7 ? t[srcEffectInd].equals("a") ? 1 : -1 : 1;
		int chg = t[srcChangeInd].startsWith("-") ? -1 : 1;
		return effect * chg;
	}
	private static int getTargetActivity(String[] t, int tgtTypeInd, int tgtEffectInd, int tgtChangeInd)
	{
		int effect = t[tgtTypeInd].endsWith("protein") && t[tgtTypeInd].length() > 7 ? (t[tgtEffectInd].equals("a") ? 1 : t[tgtEffectInd].equals("i") ? -1 : 0) : 1;
		int chg = t[tgtChangeInd].startsWith("-") ? -1 : 1;
		return effect * chg;
	}
}
