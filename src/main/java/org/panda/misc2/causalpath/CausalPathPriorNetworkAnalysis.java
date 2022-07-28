package org.panda.misc2.causalpath;

import org.panda.utility.FileUtil;
import org.panda.utility.SIFFileUtil;

import java.io.IOException;
import java.util.Set;

public class CausalPathPriorNetworkAnalysis
{
	public static void main(String[] args)
	{

	}

	public static void explorePathBetween(String priorFile, Set<String> goi, String ouFile) throws IOException
	{
		SIFFileUtil.writeSubgraph(priorFile, goi, ouFile);
	}
}
