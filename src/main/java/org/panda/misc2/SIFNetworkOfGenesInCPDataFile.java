package org.panda.misc2;

import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class SIFNetworkOfGenesInCPDataFile
{
	public static final String NETWORK_SIF_FILE = "/home/ozgunbabur/Data/PathwayCommons13.All.hgnc.txt.sif";

	static final Set<String> UNDIRECTED = new HashSet<>(Arrays.asList("interacts-with", "in-complex-with"));
	static final double LOG2 = Math.log(2);

	public static void main(String[] args) throws IOException
	{
		String base = "/home/ozgunbabur/Analyses/Kai/";

		FileUtil.processDirsRecursive(new File(base), dir ->
		{
			String data = dir + "/data.tsv";
			if (FileUtil.exists(data))
			{
				String pbName = dir + "/network-paths-between";
				String neighName = dir + "/network-neighborhood";
				Set<String> interactionTypes = new HashSet<>(Arrays.asList("controls-state-change-of",
					"controls-expression-of", "in-complex-with", "interacts-with"));

				generate(data, "ID", "Symbols", "Feature", "SignedP", 0.05,
					interactionTypes, true, pbName);
				generate(data, "ID", "Symbols", "Feature", "SignedP", 0.05,
					interactionTypes, false, neighName);
			}
		});
	}

	public static void generate(String dataFile, String idCol, String symCol, String featureCol, String signedPCol,
								double pThr, Set<String> interactionTypes, boolean pathsBetween, String outFileWOExt) throws IOException
	{
		String[] header = FileUtil.readHeader(dataFile);
		int idInd = ArrayUtil.indexOf(header, idCol);
		int pInd = ArrayUtil.indexOf(header, signedPCol);
		int symInd = ArrayUtil.indexOf(header, symCol);
		int featInd = ArrayUtil.indexOf(header, featureCol);

		Map<String, String[]> lineMap = FileUtil.linesTabbedSkip1(dataFile).collect(Collectors.toMap(t -> t[idInd], t -> t, (strings, strings2) -> strings));

		Set<String> seedSyms = new HashSet<>();
		Set<String> seedIDs = new HashSet<>();

		lineMap.forEach((id, t) ->
		{
			if (Math.abs(Double.parseDouble(t[pInd])) <= pThr)
			{
				seedSyms.add(t[symInd]);
				seedIDs.add(id);
			}
		});

		Set<String> memory = new HashSet<>();

		BufferedWriter sifWriter = FileUtil.newBufferedWriter(outFileWOExt + ".sif");

		FileUtil.linesTabbed(NETWORK_SIF_FILE).filter(t -> interactionTypes.contains(t[1])).
			filter(t -> pathsBetween ? (seedSyms.contains(t[0]) && seedSyms.contains(t[2])) :
				(seedSyms.contains(t[0]) || seedSyms.contains(t[2]))).forEach(t ->
		{
			String type = t[1];

			if (type.equals("in-complex-with") && interactionTypes.contains("interacts-with"))
			{
				type = "interacts-with";
			}

			String s = t[0] + "\t" + type + "\t" + t[2];
			if (!memory.contains(s))
			{
				FileUtil.writeln(s, sifWriter);
				memory.add(s);
				if (isUndirected(type))
				{
					s = t[2] + "\t" + type + "\t" + t[0];
					memory.add(s);
				}
			}
		});
		sifWriter.close();

		BufferedWriter fmtWriter = FileUtil.newBufferedWriter(outFileWOExt + ".format");
		fmtWriter.write("node\tall-nodes\tcolor\t255 255 255\n" +
			"node\tall-nodes\tbordercolor\t50 50 50");

		ValToColor vtc = new ValToColor(new double[]{-10, 0, 10}, new Color[]{new Color(40, 80, 255), Color.white, new Color(255, 80, 40)});
		seedIDs.forEach(id ->
		{
			String[] t = lineMap.get(id);
			double v = -Math.log(Math.abs(Double.parseDouble(t[pInd]))) / LOG2;
			if (t[pInd].startsWith("-")) v = -v;
			String color = vtc.getColorInString(v);

			if (t[featInd].equals("G"))
			{
				FileUtil.lnwrite("node\t" + t[symInd] + "\tcolor\t" + color, fmtWriter);
			}
			else
			{
				FileUtil.lnwrite("node\t" + t[symInd] + "\trppasite\t" + id + "|" + t[featInd].toLowerCase() + "|" + color + "|50 50 50|" + t[pInd], fmtWriter);
			}
		});

		fmtWriter.close();
	}

	static boolean isUndirected(String type)
	{
		return UNDIRECTED.contains(type);
	}

}
