package org.panda.misc2.teaching;

import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashMap;
import java.util.Map;

public class GradescopeSimilarityGrapher
{
	public static void main(String[] args)
	{
		generateSIFs("/home/ozgunbabur/Documents/Teaching/MachineLearning/Fall2023/ass1-sim.tsv");
	}

	public static void generateSIFs(String file)
	{
		String dir = file.substring(0, file.lastIndexOf(".")) + "/";
		FileUtil.mkdirs(dir);

		// file --> source --> target -> weight
		Map<String, Map<String, Map<String, Double>>> simMap = new HashMap<>();

		FileUtil.linesTabbed(file).forEach(t ->
		{
			String func = t[1];
			String source = t[0].split(" ")[0];
			String target = t[4].split(" ")[0];
			double weight = Double.parseDouble(t[3].substring(0, t[3].length() - 1)) / 100;

			if (!simMap.containsKey(func)) simMap.put(func, new HashMap<>());
			Map<String, Map<String, Double>> map = simMap.get(func);
			if (!map.containsKey(source)) map.put(source, new HashMap<>());
			Map<String, Double> tarMap = map.get(source);
			tarMap.put(target, weight);
		});

		for (String func : simMap.keySet())
		{
			String outSIF = dir + func + ".sif";
			String outFMT = dir + func + ".format";

			BufferedWriter sifWriter = FileUtil.newBufferedWriter(outSIF);
			BufferedWriter fmtWriter = FileUtil.newBufferedWriter(outFMT);

			FileUtil.write("graph\tgrouping\toff", fmtWriter);
			FileUtil.lnwrite("edge\tall-edges\twidth\t2", fmtWriter);

			ValToColor vtc = new ValToColor(new double[]{0.4, 1}, new Color[]{Color.white, Color.black});

			Map<String, Map<String, Double>> map = simMap.get(func);
			for (String source : map.keySet())
			{
				Map<String, Double> tarMap = map.get(source);
				for (String target : tarMap.keySet())
				{
					double weight = tarMap.get(target);

					FileUtil.writeln(source + "\tphosphorylates\t" + target, sifWriter);
					FileUtil.lnwrite("edge\t" + source + " phosphorylates " + target + "\tcolor\t" + vtc.getColorInString(weight), fmtWriter);
				}
			}

			FileUtil.closeWriter(sifWriter);
			FileUtil.closeWriter(fmtWriter);
		}
	}
}
