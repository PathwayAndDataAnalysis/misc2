package org.panda.misc2.analyses;

import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Prasidda
{
    static final String DIR = "/home/ozgunbabur/Analyses/PrasiddaPedGliomaPPM1D/";

    public static void main(String[] args) throws IOException
    {
        convert(DIR + "PPM1D-TMT6-IMAC-MedianMAD_Two-sample_mod_T_2021-09-16_n6x36506.gct", DIR + "data.tsv");
    }


    static void convert(String inFile, String outFile) throws IOException
    {
        String[] header = FileUtil.readHeader(inFile, 2);
        int sitesInd = ArrayUtil.indexOf(header, "id");
        int symbolInd = ArrayUtil.indexOf(header, "id.mapped");
        int fcInd = ArrayUtil.indexOf(header, "logFC.GSK.over.DMSO");
        int pInd = ArrayUtil.indexOf(header, "adj.P.Val.GSK.over.DMSO");

        BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
        writer.write("ID\tSymbols\tSites\tFeature\tEffect\tSignedP");

        FileUtil.linesTabbed(inFile).skip(4).forEach(t ->
        {
            String sym = t[symbolInd];

            if (sym.equals("NotFound")) return;

            String sitesSubStr = t[sitesInd].split("_")[1];
            String[] sitesArray = sitesSubStr.split("s|t|y");
            String sites = ArrayUtil.merge("|", sitesArray);

            double signedP = Double.parseDouble(t[pInd]);
            if (t[fcInd].startsWith("-")) signedP = -signedP;

            String id = sym + "-" + ArrayUtil.merge("-", sitesArray);

            FileUtil.lnwrite(id + "\t" + sym + "\t" + sites + "\tP\t\t" + signedP, writer);
        });

        writer.close();
    }
}
