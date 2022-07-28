package org.panda.misc2.causalpath;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

public class CompareResultNetworks
{
	public static void main(String[] args)
	{
		String dir = "/home/ozgunbabur/Analyses/Aslan-Thrombin-PAR/";
		compareTwo(dir + "Resting-vs-Thrombin/causative.sif", dir + "Resting-vs-Thrombin-update/causative.sif");
	}

	public static void compareTwo(String file1, String file2)
	{
		Set<String> rels1 = load(file1);
		Set<String> rels2 = load(file2);

		CollectionUtil.printNameMapping("Old", "New");
		CollectionUtil.printVennSets(50, rels1, rels2);
	}

	public static Set<String> load(String file)
	{
		return FileUtil.linesTabbed(file).filter(t -> t.length > 3)
			.map(t -> t[0] + "\t" + t[1] + "\t" + t[2]).collect(Collectors.toSet());
	}
}
