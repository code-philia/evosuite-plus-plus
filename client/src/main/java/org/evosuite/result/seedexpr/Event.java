package org.evosuite.result.seedexpr;

public abstract class Event {
	
	public static final int branchCovering = 0;
	public static final int staticPoolSampling = 1;
	public static final int staticContextPoolSampling = 2;
	public static final int dynamicPoolSampling = 3;
	public static final int randomSampling = 4;
	public static final int search = 5;
	
	private long timestamp;
	private int type;
	private String dataType;
	
	public Event(long timestamp, String dataType) {
		this.setTimestamp(timestamp);
		this.setDataType(dataType);
	}

	public long getTimestamp() {
		return timestamp;
	}

	public void setTimestamp(long timestamp) {
		this.timestamp = timestamp;
	}

	public int getType() {
		return type;
	}

	public void setType(int type) {
		this.type = type;
	}

	public String getDataType() {
		return dataType;
	}

	public void setDataType(String dataType) {
		this.dataType = dataType;
	}
}