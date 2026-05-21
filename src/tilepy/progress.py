import json
import time
from threading import Lock

_progress = {}
_lock = Lock()


def get_progress(task_id):
    with _lock:
        current_progress = _progress.get(task_id)
        # print("Initial progress fetch:", current_progress)
        # print("Current progress dict:", _progress)

    if current_progress is None:
        yield json.dumps(
            {"status": "unknown", "progress": 0, "message": "No such task_id"}
        )
    else:
        try:
            while True:
                with _lock:
                    current_progress = _progress.get(task_id)
                if current_progress is not None:
                    yield json.dumps(current_progress)
                else:
                    yield json.dumps(
                        {
                            "status": "completed",
                            "progress": 1,
                            "message": "Task completed",
                        }
                    )
                    break
                time.sleep(0.1)
        except GeneratorExit:
            print(f"Stopped monitoring progress for task_id: {task_id}")


def report(task_id, progress=None, message=None, status=None, result=None):
    if task_id is None:
        # print("Warning: report called without task_id. Progress update will be ignored.")
        return
    with _lock:
        data = _progress.get(task_id)

        if data is None:
            data = {
                "progress": 0,
                "message": "",
                "status": "in_progress",
                "time": time.time(),
            }
            print(f"Creating new progress for task_id: {task_id} at {data['time']}")

        if progress is not None:
            data["progress"] = progress
        if message is not None:
            data["message"] = message
        if status is not None:
            data["status"] = status
        print(
            f"Reporting progress for {task_id} at {time.time() - data['time']} seconds"
        )

        if status == "completed":
            print(f"Task {task_id} completed. Cleaning up progress data.")
            print(
                f"Final progress data for {task_id} in: {time.time() - data['time']} seconds"
            )
            _progress.pop(task_id, None)
        else:
            _progress[task_id] = data
        print(f"Reported progress for {task_id}: {_progress} \n")
