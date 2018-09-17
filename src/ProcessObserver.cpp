/*
 * ProcessObserver.cpp
 *
 *  Created on: Dec 6, 2013
 *      Author: mchristopher
 */

#include "ProcessObserver.h"

ProcessObserver::ProcessObserver(){

	type = OptIteration;
	count = 0;

}

void ProcessObserver::Execute(itk::Object *caller, const itk::EventObject & event){
	Execute( (const itk::Object *)caller, event);
}

void ProcessObserver::Execute(const itk::Object * object, const itk::EventObject & event){
	if(this->type == OptIteration && itk::IterationEvent().CheckEvent(&event)){

		const itk::SingleValuedNonLinearOptimizer *optimizer =
				dynamic_cast<const itk::SingleValuedNonLinearOptimizer *>(object);

		std::cout << count << " = ";
		std::cout << optimizer->GetValue(optimizer->GetCurrentPosition());
		std::cout << " (" << optimizer->GetCurrentPosition() << ")";
		std::cout << std::endl;

//		std::cout << count << " = ";
//		std::cout << opt->GetValue(opt->GetCurrentPosition());
//		std::cout << " (" << opt->GetCurrentPosition() << ")";
//		std::cout << std::endl;

		std::cout.flush();
		++count;
	}
}

