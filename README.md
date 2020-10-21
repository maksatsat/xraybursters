# X-ray bursters

## Генератор событий  

Принимаемые параметры - `m, intervals`

`m` - число событий в единицу времени

`intervals` - список интервалов по времени

```
from event_generator import Events
a = Events(1/3, [[100, 300], [400,500]])
```

![example](example.png)