# skele/commands/base.py
"""The base command."""


class Cmd:
    """A base command."""

    def __init__(self):
        pass

    def run(self):
        raise NotImplementedError('You must implement the run() method yourself!')

