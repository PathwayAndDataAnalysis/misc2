package org.panda.misc2;

import org.panda.causalpath.network.GraphWriter;
import org.panda.utility.FileUtil;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Set;

public class CopySIFFilesRecursively
{
	public static void main(String[] args) throws IOException
	{
		String inBase = "/home/ozgunbabur/Analyses/CPTAC-PanCan/";
		String inDir = inBase + "mutational-signatures/";
		String outBase = inBase + "sif-files/";

		generate(inBase, inDir, outBase);
	}

	public static void generate(String inBase, String inDir, String outBase) throws IOException
	{
		for (File file : new File(inDir).listFiles())
		{
			if (file.getPath().endsWith(".sif"))
			{
				File format = new File(file.getPath().substring(0, file.getPath().lastIndexOf(".")) + ".format");

				if (format.exists())
				{
					String outDir = inDir.replace(inBase, outBase);
					Files.createDirectories(Paths.get(outDir));

					FileUtil.copyFile(file.getPath(), file.getPath().replace(inBase, outBase));
					FileUtil.copyFile(format.getPath(), format.getPath().replace(inBase, outBase));
				}
			}
		}

		for (File sub : new File(inDir).listFiles())
		{
			if (sub.isDirectory())
			{
				generate(inBase, sub.getPath(), outBase);
			}
		}
	}
}
