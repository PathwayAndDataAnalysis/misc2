package org.panda.misc2;

import org.panda.utility.FileUtil;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class ScPerturbStats
{
	public static final String DIR = "/home/ozgunbabur/Data/scPerturb/summary/";
	public static void main(String[] args)
	{
		printCommonPerturb();
	}

	public static void printCommonPerturb()
	{
		Map<String, Map<String, Integer>> common = new HashMap<>();
		for (File file : new File(DIR).listFiles())
		{
			if (file.getName().contains("perturbation_counts"))
			{
				String[] s = file.getName().split("[/_.]");
				String study = s[0] + "_" + s[1];

				System.out.println("study = " + study);

				FileUtil.linesTabbedSkip1(file.getPath()).forEach(t ->
				{
					if (!common.containsKey(t[0])) common.put(t[0], new HashMap<>());
					common.get(t[0]).put(study, Integer.parseInt(t[1]));
				});
			}
		}

		Map<String, Set<String>> priors = loadPriors();

		common.forEach((target, map) ->
		{
			if (priors.keySet().contains(target) && priors.get(target).size() >= 10 && map.size() > 1)
			{
				System.out.println("\ntarget = " + target + "\tnum targs = " + priors.get(target).size());
				map.forEach((study, cnt) -> System.out.println(study + "\t" + cnt));
			}
		});
	}

	public static Map<String, Set<String>> loadPriors()
	{
		Map<String, Set<String>> priors = new HashMap<>();
		FileUtil.linesTabbed("/home/ozgunbabur/Data/causal-priors.txt")
			.filter(t -> t[1].equals("upregulates-expression") || t[1].equals("downregulates-expression")).forEach(t ->
		{
			if (!priors.containsKey(t[0])) priors.put(t[0], new HashSet<>());
			priors.get(t[0]).add(t[2]);
		});
		return priors;
	}
}
