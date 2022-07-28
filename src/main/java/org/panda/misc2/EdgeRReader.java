package org.panda.misc2;

public class EdgeRReader
{
	public static String extractGeneSymbolFromDescription(String[] t, int symInd)
	{
		String s = t[symInd];
		int start = s.indexOf(" GN=");
		if (start >= 0)
		{
			int end = s.indexOf(" ", start + 4);

			if (end > 0) return s.substring(start + 4, end);
			else return s.substring(start + 4);
		}

		return null;
	}

	public static double readSignedP(String[] t, int fcInd, int pInd)
	{
		if (!t[fcInd].isEmpty() && !t[pInd].isEmpty())
		{
			double fc = Double.parseDouble(t[fcInd]);
			double p = Double.parseDouble(t[pInd]);

			if (Double.isNaN(p)) return Double.NaN;
			else if (p == 0)
			{
				p = 1e-20;
			}

			return Math.signum(fc) * p;
		}

		return Double.NaN;
	}
}
