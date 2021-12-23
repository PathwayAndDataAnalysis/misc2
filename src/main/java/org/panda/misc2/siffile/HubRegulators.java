package org.panda.misc2.siffile;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.SIFFileUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class HubRegulators
{
    public static void writeHubRegulatorsRecursive(String srcDir, String trgDir, String sifName, String outName) throws IOException
    {
        FileUtil.processDirsRecursive(new File(srcDir), dir ->
        {
            String inFile = dir.getPath() + File.separator + sifName;
            if (FileUtil.exists(inFile))
            {
                String outFile = dir.getPath().replaceFirst(srcDir, trgDir) + File.separator + outName;
                writeHubRegulators(inFile, outFile);
            }
        });
    }
    public static void writeHubRegulators(String sifFile, String outFile) throws IOException
    {
        FileUtil.mkdirsOfFilePath(outFile);
        BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
        writer.write("Regulator\tTarg size\tTargets");
        Map<String, Set<String>> map = SIFFileUtil.convertToMap(sifFile);
        map.keySet().stream().sorted((s1, s2) -> Integer.compare(map.get(s2).size(), map.get(s1).size())).forEach(s ->
        {
            List<String> dwns = new ArrayList<>(map.get(s));
            dwns.sort(String::compareTo);

            FileUtil.lnwrite(s + "\t" + dwns.size() + "\t" + CollectionUtil.merge(dwns, " "), writer);
        });
        writer.close();
    }
}
