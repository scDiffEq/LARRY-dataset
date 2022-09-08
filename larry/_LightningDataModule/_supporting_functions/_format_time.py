
def _format_time(self, task, train_key, test_key, time_key):

    """format time"""

    self._task = task
    self._train_key = train_key
    self._test_key = test_key
    self._time_key = time_key

    time_dict = {
        "timepoint_recovery": {self._train_key: [2, 6], self._test_key: [2, 4]},
        "fate_prediction": {self._train_key: [2, 4, 6], self._test_key: [2, 4, 6]},
    }

    self._train_time = time_dict[self._task][self._train_key]
    self._test_time = time_dict[self._task][self._test_key]