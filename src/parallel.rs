pub trait Parallel: Sync {
    fn join<A: Send, B: Send>(
        &self,
        a: impl FnOnce() -> A + Send,
        b: impl FnOnce() -> B + Send,
    ) -> (A, B);

    fn for_each<T: Send>(&self, items: &mut [T], job: impl Fn(usize, &mut T) + Sync);
}

pub trait ForEach<T>: Sync {
    fn each(&self, items: &mut [T], job: &(dyn Fn(usize, &mut T) + Sync));
}

impl<T: Send, P: Parallel> ForEach<T> for P {
    fn each(&self, items: &mut [T], job: &(dyn Fn(usize, &mut T) + Sync)) {
        self.for_each(items, |index, item| job(index, item));
    }
}

pub struct Sequential;

impl Parallel for Sequential {
    fn join<A: Send, B: Send>(
        &self,
        a: impl FnOnce() -> A + Send,
        b: impl FnOnce() -> B + Send,
    ) -> (A, B) {
        (a(), b())
    }

    fn for_each<T: Send>(&self, items: &mut [T], job: impl Fn(usize, &mut T) + Sync) {
        for (index, item) in items.iter_mut().enumerate() {
            job(index, item);
        }
    }
}
