package org.panda.misc2.causalpath;

import org.panda.resource.ChEBI;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;

public class ChEBINameSwitcher
{
	public static final String ID_PREFIX = "CHEBI:";

	public static void main(String[] args)
	{
		String dir = "/home/ozgunbabur/Analyses/Platelet-Blood-paper/ghatge-metabolome/cond2-relax1aa/";
		replaceCheBIIDsWithNames(dir + "BTK-downstream", dir + "BTK-downstream-withchemnames");
	}

	public static void replaceCheBIIDsWithNames(String inputSIFWOExt, String outSIFWOExt)
	{
		replaceChEBIIDsOnEachLine(inputSIFWOExt + ".sif", outSIFWOExt + ".sif");
		if (FileUtil.exists(inputSIFWOExt + ".format"))
		{
			replaceChEBIIDsOnEachLine(inputSIFWOExt + ".format", outSIFWOExt + ".format");
		}
		if (FileUtil.exists(inputSIFWOExt + ".layout"))
		{
			replaceChEBIIDsOnEachLine(inputSIFWOExt + ".layout", outSIFWOExt + ".layout");
		}
	}

	private static void replaceChEBIIDsOnEachLine(String inFile, String outFile)
	{
		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);

		FileUtil.lines(inFile).forEach(l ->
		{
			int[] loc = getLocationOfAChEBIID(l);
			while (loc != null)
			{
				String id = l.substring(loc[0], loc[1]);
				String name = ChEBI.get().getName(id);
				if (name != null)
				{
					l = l.substring(0, loc[0]) + name + l.substring(loc[1]);
					loc = getLocationOfAChEBIID(l);
				}
				else
				{
					System.err.println("No name for: " + id);
					break;
				}
			}
			FileUtil.writeln(l, writer);
		});

		FileUtil.closeWriter(writer);
	}

	private static int[] getLocationOfAChEBIID(String line)
	{
		int start = line.indexOf(ID_PREFIX);
		if (start >= 0)
		{
			int end = getIndexOfNonDigit(line, start + 6);
			return new int[]{start, end};
		}
		else return null;
	}

	private static int getIndexOfNonDigit(String line, int start)
	{
		int i = start;
		while (line.length() > i && Character.isDigit(line.charAt(i))) i++;
		return i;
	}
}
