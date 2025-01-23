package org.panda.misc2;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Progress;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;

import java.io.BufferedWriter;
import java.util.*;

public class TestingMeanRank
{
	static Random rand = new Random();

	public static void main(String[] args)
	{
//		int n = 100;
//		int k = 20;
//		int size = 1000000;
//		double[] distribution = getDistribution(n, k, size);
//		double stdev = Summary.stdev(distribution);
//		System.out.println("stdev = " + stdev);

		double p = (1 + Erf.erf(1 / Math.sqrt(2))) / 2;
		System.out.println("p = " + p);


//		testMany();
	}

	public static void testMany()
	{
		int n = 100;
		int k = 10;
		int size = 100000000;
//		for (int k = 5; k < 50; k++)

//		for (int n = 10; n <= 10000; n++)
//		{
//			List<Double> distribution = getDistribution(n, k, size);
//			double p = getPForThr(distribution, 0.3);
//			System.out.println(n + "\t" + p);
//		}

		Histogram cumDist = new Histogram(0.01);
		cumDist.setBorderAtZero(true);
		Histogram cumCtrl = new Histogram(0.01);
		cumCtrl.setBorderAtZero(true);

		Histogram histKS = new Histogram(0.1);
		Histogram histPDiff = new Histogram(0.0001);
		histKS.setBorderAtZero(true);

		KolmogorovSmirnovTest ks = new KolmogorovSmirnovTest();

//		Progress prg = new Progress(1000, "Calculating distributions");
		for (int i = 0; i < 2; i++)
		{
//			prg.tick();
			double[] distribution = getDistribution(n, k, size);


			cumDist.countAll(distribution);
//			double mean = Summary.mean(distribution);
//			System.out.println("mean = " + mean);
			double stdev = Summary.stdev(distribution);
			System.out.println(n + "\t" + stdev);

			double[] ctrlDist = new double[size];
			for (int j = 0; j < size; j++)
			{
				ctrlDist[j] = (rand.nextGaussian() * stdev) + 0.5;
				cumCtrl.count(ctrlDist[j]);
			}
//			stdev = Summary.stdev(ctrlDist);
//			System.out.println("stdev = " + stdev);

//			double ksPval = ks.kolmogorovSmirnovTest(distribution, ctrlDist);
//			histKS.count(ksPval);

//			double pTest = getPForThr(distribution, 0.35);
//			double pCtrl = getPForThr(ctrlDist, 0.35);
//			double dif = pTest - pCtrl;
//			histPDiff.count(dif);
		}

		cumDist.printTogether(cumCtrl);
//		histKS.print();

		System.out.println(" -------------------- ");

//		histPDiff.print();

//		double mean = Summary.meanOfDoubles(distribution);
//		double stdev = Summary.stdevOfDoublesInList(distribution);
//		System.out.println("mean = " + mean + "\tstdev = " + stdev);


//		Histogram h = new Histogram(0.01);
//		h.setBorderAtZero(true);
//		distribution.forEach(h::count);
//		h.print();
	}

	public static double[] getDistribution(int n, int k, int size)
	{
		double[] distr = new double[size];

		for (int i = 0; i < size; i++)
		{
			Set<Double> set = new HashSet<>();

//			for (int j = 0; j < k; j++)
			while (set.size() < k)
			{
				double rank = (rand.nextInt(n) + 0.5) / n;
				set.add(rank);
			}

			double mr = Summary.sum(set) / k;

//			if (mr > 0.5) mr = 1 - mr;

			distr[i] = mr;
		}
		return distr;
	}

	static double getPForThr(double[] distr, double thr)
	{
		int cnt = 0;
		for (double v : distr)
		{
			if (v <= thr) cnt++;
		}
		return cnt / (double) distr.length;
	}
}
