import java.awt.*;
import java.awt.print.*;

/**
 *PrintableComponent.java
 *
 *Created on 26.02.2001
 *
 *Last modified on 26.02.2001
 */


public class PrintableComponent implements Printable
{

Component c;

public PrintableComponent(Component c)
	{
	this.c = c;
	}

public void print() throws PrinterException
	{
	PrinterJob job = PrinterJob.getPrinterJob();
	PageFormat format = job.pageDialog(job.defaultPage());
	job.setPrintable(this,format);
	if(job.printDialog())
		job.print();
	}

public int print(Graphics g, PageFormat format, int pagenum)
	{
	if(pagenum>0)
		{
		return Printable.NO_SUCH_PAGE;
		}
	Graphics2D g2 = (Graphics2D)g;
	g2.translate(format.getImageableX(),format.getImageableY());
	c.paint(g2);
	return Printable.PAGE_EXISTS;
	}

}