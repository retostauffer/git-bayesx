import javax.swing.SwingUtilities;

/**
 *ParseThread.java
 *
 *Created on 21.03.2001
 *
 *Last modified on 21.03.2001
 */


public class ParseThread
{

private Object value;
private Thread thread;
private BayesX b;
private String input;

private static class ThreadVar
	{
	private Thread thread;
	ThreadVar(Thread t)
		{
		thread = t;
		}
	synchronized Thread get()
		{
		return thread;
		}
	synchronized void clear()
		{
		thread = null;
		}
	}

private ThreadVar threadVar;


protected synchronized Object getValue()
	{
	return value; 
	}

private synchronized void setValue(Object x)
	{
	value = x; 
	}

public Object construct()
	{
	Boolean st = new Boolean(b.parse(input));
	return st;
	}


public void finished() 
	{
	boolean stop = ((Boolean)this.get()).booleanValue();
	b.setProcessrunning(false);
	b.finished();
	if(stop)
		{
		b.fileCommand(1);
		}
	}

public void interrupt()
	{
	Thread t = threadVar.get();
	if (t != null)
		{
		t.interrupt();
		}
	threadVar.clear();
	}


public Object get()
	{
	while (true)
		{
		Thread t = threadVar.get();
		if (t == null)
			{
			return getValue();
			}
		try
			{
			t.join();
			}
		catch (InterruptedException e)
			{
			Thread.currentThread().interrupt();
			return null;
			}
		}
	}


public void setCommand(String str)
	{
	input = str;
	}

public void setPriority(int p)
	{
	Thread t = threadVar.get();
	if (t!= null)
		{
		t.setPriority(p);
		}
	}

public int getPriority()
	{
	Thread t = threadVar.get();
	if (t!= null)
		{
		return t.getPriority();
		}
	else
		{
		return 0;
		}
	}

/*public void suspend()
	{
	Thread t = threadVar.get();
	if (t != null)
		{
		t.suspend();
		}
	}

public void resume()
	{
	Thread t = threadVar.get();
	if (t != null)
		{
		t.resume();
		}
	}*/

public ParseThread(BayesX b)
	{
	this.b=b;
	final Runnable doFinished = new Runnable()
		{
		public void run()
			{
			finished();
			}
		};
	Runnable doConstruct = new Runnable()
		{ 
		public void run()
			{
			try
				{
				setValue(construct());
				}
			finally
				{
				threadVar.clear();
				}
			SwingUtilities.invokeLater(doFinished);
			}
		};
	Thread t = new Thread(doConstruct);
	threadVar = new ThreadVar(t);
//	setPriority(Thread.NORM_PRIORITY);
	}

public void start()
	{
	Thread t = threadVar.get();
	if (t != null)
		{
		t.start();
		}
	}

}
