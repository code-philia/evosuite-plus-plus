package com.test;

public class Example {
	public boolean test(int a, int b){
		if(Util.isOk(a, b)){
			if(Util.isCornerCase(a, b)){
				return true;
			}	
		}
		
		return false;
	}
}