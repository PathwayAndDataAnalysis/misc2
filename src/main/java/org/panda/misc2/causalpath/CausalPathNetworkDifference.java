package org.panda.misc2.causalpath;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class CausalPathNetworkDifference
{
	public static void produceDiff(String inputSIF, String[] excludeSIFs, String outSIF) throws IOException
	{
		Map<String, String> inputMap = new HashMap<>();
		FileUtil.lines(inputSIF).forEach(l ->
		{
			String[] t = l.split("\t");
			if (t.length < 3) return;
			String key = t[0] + "\t" + t[1] + "\t" + t[2];
			inputMap.put(key, l);
		});

		Set<String> exclude = new HashSet<>();
		for (String excludeSIF : excludeSIFs)
		{
			exclude.addAll(loadSIFKeys(excludeSIF));
		}

		BufferedWriter sifWriter = FileUtil.newBufferedWriter(outSIF);
		inputMap.keySet().stream().filter(k -> !exclude.contains(k)).forEach(k ->
			FileUtil.writeln(inputMap.get(k), sifWriter));

		sifWriter.close();

		String inputFmt = toFmt(inputSIF);
		if (FileUtil.exists(inputFmt)) FileUtil.copyFile(inputFmt, toFmt(outSIF));

		String inputLyt = toLyt(inputSIF);
		if (FileUtil.exists(inputLyt)) FileUtil.copyFile(inputLyt, toLyt(outSIF));
	}

	public static Set<String> loadSIFKeys(String file)
	{
		return FileUtil.linesTabbed(file).filter(t -> t.length > 2).map(t -> t[0] + "\t" + t[1] + "\t" + t[2]).collect(Collectors.toSet());
	}

	static String toFmt(String sifFilename)
	{
		return sifFilename.substring(0, sifFilename.lastIndexOf(".")) + ".format";
	}
	static String toLyt(String sifFilename)
	{
		return sifFilename.substring(0, sifFilename.lastIndexOf(".")) + ".layout";
	}

	public static void main(String[] args) throws IOException
	{
		String dir = "/home/ozgunbabur/Analyses/Aslan-Thrombin-PAR/Regular-CP/relax-2aa/";
		String par1 = dir + "Resting-vs-PAR1/causative.sif";
		String par4 = dir + "Resting-vs-PAR4/causative.sif";
		String thrombin = dir + "Resting-vs-Thrombin/causative.sif";

		produceDiff(par1, new String[]{par4, thrombin}, dir + "diff-graphs/PAR1-only.sif");
		produceDiff(par4, new String[]{par1, thrombin}, dir + "diff-graphs/PAR4-only.sif");
		produceDiff(thrombin, new String[]{par1, par4}, dir + "diff-graphs/Thrombin-only.sif");

		produceDiff(par1, new String[]{par4}, dir + "diff-graphs/PAR1-diff-PAR4.sif");
		produceDiff(par4, new String[]{par1}, dir + "diff-graphs/PAR4-diff-PAR1.sif");

	}
}
