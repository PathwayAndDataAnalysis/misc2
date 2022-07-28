package org.panda.misc2.siffile;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

public class DownstreamPathExtractor
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/home/ozgunbabur/Analyses/Platelet-Blood-paper/merged/";
		extract(dir + "merged-relax1aa.sif", dir + "GP6-downstream.sif", "GP6");
	}

	public static void extract(String inFile, String outFile, String... sources) throws IOException
	{
		Set<String> lines = FileUtil.lines(inFile).filter(l -> l.split("\t").length > 2).collect(Collectors.toSet());

		Set<String> subsetLines = new HashSet<>();

		Set<String> protsToFollow = new HashSet<>();
		Set<String> nextProtsToFollow = new HashSet<>(Arrays.asList(sources));
		Set<String> visited = new HashSet<>();

		while (!nextProtsToFollow.isEmpty())
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
}
