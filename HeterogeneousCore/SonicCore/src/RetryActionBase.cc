#include "HeterogeneousCore/SonicCore/interface/RetryActionBase.h"

// Constructor implementation
RetryActionBase::RetryActionBase(const edm::ParameterSet& conf, SonicClientBase* client)
    : client_(client), shouldRetry_(true) {}

void RetryActionBase::eval() {
  if (client_) {
    client_->evaluate();
  } else {
    edm::LogError("RetryActionBase") << "Client pointer is null, cannot evaluate.";
  }
}

void RetryActionBase::finish(bool success) {
  if (client_) {
    client_->finish(success);
  } else {
    edm::LogError("RetryActionBase") << "Client pointer is null, cannot call finish.";
  }
}

EDM_REGISTER_PLUGINFACTORY(RetryActionFactory, "RetryActionFactory");
