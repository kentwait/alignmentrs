from datetime import datetime
import json


__all__ = ['History', 'Record']


class History:
    def __init__(self, *args, **kwargs):
        self.records = []
        # self.store_state_history = store_state_history

    def add(self, operation, after, args=None, kwargs=None,
            # before_state=None, after_state=None
            ):
        record = Record(operation, after, args=args, kwargs=kwargs)
        self.records.append(record)
    
    def to_json(self, path=None):
        return '[{}]'.format(
            ', '.join([i.to_json() for i in self.records])
        )

    def to_markdown(self, path=None):
        return '\n'.join([i.to_markdown() for i in self.records])

    def to_string(self, str_formatter=None, strftime='%m/%d/%Y %H:%M:%S',
                  follow_threshold=False):
        return '\n'.join([
            i.to_string(str_formatter=str_formatter, strftime=strftime,
                        follow_threshold=follow_threshold)
            for i in self.records
        ])

    def to_toml(self, path=None):
        return '\n'.join([i.to_toml() for i in self.records])

    def to_csv(self, path=None, sep='\t'):
        return '\n'.join([i.to_csv(sep=sep) for i in self.records])

    def __str__(self):
        return '\n'.join([str(item) for item in self.records])

    def __repr__(self):
        return repr(self.records)

    def __getitem__(self, k):
        return self.record[k]

    def __setitem__(self, k, v):
        self.records[k] = v


class Record:
    def __init__(self, operation: str, state: dict,
                 args: list=None, kwargs: dict=None,
                 str_formatter=None, module=None):
        self.operation = operation
        self.optype = 'operation'
        if module:
            self.operation = '.'.join([module, self.operation])
        self.args = [repr(arg).replace('\n', ' ') for arg in args] \
            if args else []
        self.kwargs = {key: repr(kwarg).replace('\n', ' ') 
                       for key, kwarg in kwargs.items()} \
            if kwargs else dict()
        self.state = state
        self.datetime = datetime.now()
        self._threshold_args = 3
        self._threshold_kwargs = 2
        self.str_formatter = str_formatter

    def to_json(self, path=None):
        return json.dumps(
            {self.optype: 
                {
                    'datetime': self.datetime.strftime('%m/%d/%Y %H:%M:%S'),
                    'operation': self.operation,
                    'args': self.args,
                    'kwargs': self.kwargs,
                    'state': self.state,
                }
            }
        )

    def to_markdown(self, path=None):
        params = ''
        if len(self.args) > 0:
            params += ', '.join(self.args)
        if len(self.kwargs) > 0:
            if len(self.args) > 0:
                params += ', '
            params += ', '.join(
                ['{}={}'.format(k, v) for k, v in self.kwargs.items()]
            )
        return ('## {optype} {op}\n'
                '  * date and time - {dt}\n'
                '  * statement - `{op}({params})`\n'
                '  * state - `{state}`\n'.format(
                   dt=self.datetime.strftime('%m/%d/%Y %H:%M:%S'),
                   op=self.operation,
                   optype=self.optype,
                   params=params,
                   state=str(self.state),
                ))

    def to_string(self, str_formatter=None, strftime='%m/%d/%Y %H:%M:%S',
                  follow_threshold=False):
        if str_formatter:
            return str_formatter(
                operation=self.operation,
                optype=self.optype,
                args=self.args,
                kwargs=self.kwargs,
                datetime=self.datetime.strftime(strftime),
            )
        if self.str_formatter:
            return self.str_formatter(
                operation=self.operation,
                optype=self.optype,
                args=self.args,
                kwargs=self.kwargs,
                datetime=self.datetime,
            )
        if follow_threshold:
            args = [str(v) for i, v in enumerate(self.args)
                    if i < self._threshold_args]
            if len(self.args) > self._threshold_args:
                args += ['...']
            
            kwargs = ['{k}={v}'.format(k=kv[0], v=kv[1]) 
                      for i, kv in enumerate(self.kwargs.items())
                      if i < self._threshold_args]
            if len(self.kwargs) > self._threshold_kwargs:
                kwargs += ['...']
        else:
            args = [str(v) for v in self.args]
            kwargs = ['{}={}'.format(k, v) for k, v in self.kwargs.items()]
        
        params = ''
        if len(args) > 0:
            params += ', '.join(args)
        if len(kwargs) > 0:
            if len(self.args) > 0:
                params += ', '
            params += ', '.join(kwargs)
        return '{optype} {dt} {op}({params})'.format(
            dt=self.datetime.strftime('%m/%d/%Y %H:%M:%S'),
            optype=self.optype,
            op=self.operation,
            params=params,
        )

    def to_toml(self, path=None):
        return ('[{optype}]\n'
                'datetime = {dt}\n'
                'operation = \'{op}\'\n'
                'args = \'{args}\'\n'
                'kwargs = \'{kwargs}\'\n'
                'state = {state}\n'.format(
                   dt=self.datetime.strftime('%m/%d/%Y %H:%M:%S'),
                   op=self.operation,
                   optype=self.optype,
                   args=self.args,
                   kwargs=str(self.kwargs),
                   state=str(self.state),
               ))

    def to_csv(self, path=None, sep='\t'):
        return sep.join([
            self.optype,
            self.datetime.strftime('%m/%d/%Y %H:%M:%S'),
            self.operation,
            str(self.args),
            str(self.kwargs),
            str(self.state),
        ])

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return '<{clsname} datetime={dt} statement={op}({params})>'.format(
            clsname=self.__class__.__name__,
            dt=self.datetime.strftime('%m/%d/%YT%H:%M:%S'),
            op=self.operation,
            params='...' if len(self.args) > 0 or len(self.kwargs) > 0 else ''
        )
