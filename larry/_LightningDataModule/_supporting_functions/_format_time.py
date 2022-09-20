
def _format_time(self, task, train_key, test_key, time_key):

    """format time"""

    time_dict = {
        "timepoint_recovery": {train_key: [2, 6], test_key: [2, 4]},
        "fate_prediction": {train_key: [2, 4, 6], test_key: [2, 4, 6]},
    }

    self.train_time = time_dict[task][train_key]
    self.test_time  = time_dict[task][test_key]