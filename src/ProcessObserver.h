/*
 * ProcessObserver.h
 *
 *  Created on: Dec 6, 2013
 *      Author: mchristopher
 */

#include "itkCommand.h"
#include "itkSingleValuedNonLinearOptimizer.h"

#include "itkConjugateGradientOptimizer.h"

#ifndef PROCESSOBSERVER_H_
#define PROCESSOBSERVER_H_

enum EventType{
	Start = 0,
	End,
	OptIteration,
};

class ProcessObserver : public itk::Command {

private:
	/** Type of event monitored by this observer */
	EventType type;
	int count;

	itk::SingleValuedNonLinearOptimizer::Pointer opt;

public:

	ProcessObserver();
	~ProcessObserver(){};

	itkNewMacro( ProcessObserver );

	void Execute(itk::Object *caller, const itk::EventObject & event);
	void Execute(const itk::Object * object, const itk::EventObject & event);

	inline void setOptimizer(itk::SingleValuedNonLinearOptimizer::Pointer o){
		this->opt = o;
	}

};

#endif /* PROCESSOBSERVER_H_ */
