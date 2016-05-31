/***************************************************************************
 *  Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#include "gmock/gmock.h"

#include "lib/definitions.h"
#include "lib/datastructure/StateVector.h"

using::testing::Eq;
using::testing::DoubleEq;
using::testing::Test;

namespace datastructure {
class AStateVector : public Test {
 public:
  AStateVector() :
    state_vector(20) { }

  StateVector<uint16_t, 20> state_vector;
};


TEST_F(AStateVector, ReturnInvalidStateIfNoEntryIsSet) {
  for(int i = 0; i < 20; ++i) {
     ASSERT_FALSE(state_vector.isEntryValid(i));
     ASSERT_FALSE(state_vector[i]); 
  }
}

TEST_F(AStateVector, StoresCorrectStates) {
  for(int i = 0; i < 20; ++i) {
     state_vector.setState(i,i+1);
     ASSERT_THAT(state_vector[i],Eq(i+1));
  }
}

TEST_F(AStateVector, InvalidesAllEntriesAfterReset) {
  for(int i = 0; i < 20; ++i) {
     state_vector.setState(i,i+1);
  }
  state_vector.reset();
  for(int i = 0; i < 20; ++i) {
     ASSERT_FALSE(state_vector.isEntryValid(i));
     ASSERT_FALSE(state_vector[i]); 
  }
}

TEST_F(AStateVector, WorksEqualAfterReset) {
  state_vector.reset();
  for(int i = 0; i < 20; ++i) {
     state_vector.setState(i,i+1);
     ASSERT_TRUE(state_vector.isEntryValid(i));
     ASSERT_THAT(state_vector[i],Eq(i+1));
  }
}

TEST_F(AStateVector, WorksEqualAfterMultipleResetOperations) {
  for(int j = 0; j < 5; ++j) {
    state_vector.reset();
    for(int i = 0; i < 20; ++i) {
     ASSERT_FALSE(state_vector.isEntryValid(i));
     ASSERT_FALSE(state_vector[i]); 
    }
    for(int i = 0; i < 20; ++i) {
      state_vector.setState(i,i+1);
      ASSERT_TRUE(state_vector.isEntryValid(i));
      ASSERT_THAT(state_vector[i],Eq(i+1));
    }
  }
}


}  // namespace datastructure
