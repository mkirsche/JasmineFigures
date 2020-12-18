import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.Scanner;

public class AddStrandBiasHeader {
public static void main(String[] args) throws Exception
{
	String fn = args[0], ofn = args[1];
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	PrintWriter out = new PrintWriter(new File(ofn));
	boolean headerDone = false;
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0)
		{
			continue;
		}
		if(!line.startsWith("##"))
		{
			if(!headerDone)
			{
				out.println("##FILTER=<ID=STRANDBIAS,Description=\"Strand bias detected by Sniffles.\">");
				out.println("##FILTER=<ID=PASS,Description=\"Pass detected by Sniffles.\">");
			}
			headerDone = true;
		}
		if(line.startsWith("##INFO=<ID=INTRASAMPLE_IDLIST"))
		{
			out.println("##INFO=<ID=INTRASAMPLE_IDLIST,Number=1,Type=String,Description=\"The IDs which were merged in the most recent round of merging\">");
		}
		else
		{
			out.println(line.replaceAll(";SUPP_VE=", ";SUPP_VEC="));
		}
	}
	out.close();
}
}
